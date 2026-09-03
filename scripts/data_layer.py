# -*- coding: utf-8 -*-
"""Canonical data layer for mechanism-oriented ADNI exploration."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C
from trajectory import per_patient_change

GROUP_LABELS = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}
ALIGNMENT_WINDOWS = {"90d": 90, "180d": 180, "365d": 365}


@dataclass(frozen=True)
class AlignmentSpec:
    name: str
    days: int


def _read_csv(key: str) -> pd.DataFrame:
    path = C.RAW_DIR / C.FILES[key]
    if not path.exists():
        raise FileNotFoundError(f"缺少原始文件: {path}")
    return pd.read_csv(path, low_memory=False)


def clean_numeric(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return values.mask(values <= 0)


def earliest_row(df: pd.DataFrame, date_col: str, value_cols: list[str] | None = None,
                 rid_col: str = "RID") -> pd.DataFrame:
    """Select one source row per RID; never use column-wise groupby().first()."""
    value_cols = value_cols or []
    if rid_col not in df.columns or date_col not in df.columns:
        raise KeyError(f"需要列 {rid_col}/{date_col}")
    source = df.copy()
    source[date_col] = pd.to_datetime(source[date_col], errors="coerce")
    source = source.dropna(subset=[rid_col, date_col]).sort_values([rid_col, date_col])
    return source.drop_duplicates(rid_col, keep="first")


def _baseline_table(key: str, value_cols: list[str]) -> pd.DataFrame:
    table = _read_csv(key)
    return earliest_row(table, C.FIELDS["date"][key], value_cols)


def build_baseline() -> pd.DataFrame:
    """Build same-row baseline data while preserving missingness."""
    demo = _baseline_table("ptdemog", [C.FIELDS["pt_dob_yy"], C.FIELDS["gender"], C.FIELDS["educat"]])
    date_col = C.FIELDS["date"]["ptdemog"]
    visit_date = pd.to_datetime(demo[date_col], errors="coerce")
    # PTDOBYY 是出生年份 float（如 1931.0）：不能喂给 pd.to_datetime（会被当作纳秒时间戳），年份数值直接相减。
    dob_yy = pd.to_numeric(demo[C.FIELDS["pt_dob_yy"]], errors="coerce")
    demo["AGE"] = visit_date.dt.year - dob_yy
    out = demo[["RID", "AGE", C.FIELDS["gender"], C.FIELDS["educat"]]].rename(
        columns={C.FIELDS["gender"]: "GENDER", C.FIELDS["educat"]: "EDUCAT"})
    out["GENDER"] = C.normalize_gender(out["GENDER"])  # ADNIMERGE2 字符串 → 数字编码

    for key, source_col, output_col in [
        ("cdr", C.FIELDS["cdrsb"], "CDRSB_bl"),
        ("mmse", C.FIELDS["mmse"], "MMSE_bl"),
        ("adas", C.FIELDS["adas13"], "ADAS13_bl"),
        ("dxsum", C.FIELDS["dx"], "DX_bl"),
    ]:
        table = _baseline_table(key, [source_col])
        if key == "dxsum":
            values = C.normalize_dx(table[source_col])  # ADNIMERGE2 'CN'/'MCI'/'Dementia' → 1/2/3
        else:
            values = pd.to_numeric(table[source_col], errors="coerce")
        part = pd.DataFrame({"RID": table["RID"], output_col: values})
        out = out.merge(part, on="RID", how="left")

    plasma = _baseline_table("plasma", [C.FIELDS["p_tau217"], C.FIELDS["abeta42_ab40"]]).copy()
    plasma["PTAU217_bl"] = clean_numeric(plasma[C.FIELDS["p_tau217"]])
    plasma["AB_RATIO_bl"] = clean_numeric(plasma[C.FIELDS["abeta42_ab40"]])
    out = out.merge(plasma[["RID", "PTAU217_bl", "AB_RATIO_bl"]], on="RID", how="left")

    pet = _read_csv("pet")
    pet[C.FIELDS["pet_status"]] = pd.to_numeric(pet[C.FIELDS["pet_status"]], errors="coerce")
    pet = pet.dropna(subset=[C.FIELDS["pet_status"]])
    pet = earliest_row(pet, C.FIELDS["date"]["pet"], [C.FIELDS["pet_status"]])
    pet_cols = ["RID", C.FIELDS["pet_status"], C.FIELDS["centiloids"], C.FIELDS["tracer"]]
    out = out.merge(pet[pet_cols].rename(columns={
        C.FIELDS["pet_status"]: "PET_STATUS",
        C.FIELDS["centiloids"]: "CENTILOIDS_bl",
        C.FIELDS["tracer"]: "PET_TRACER",
    }), on="RID", how="left")
    out["PET_STATUS"] = pd.to_numeric(out["PET_STATUS"], errors="coerce")
    return out


def compute_plasma_status(df: pd.DataFrame, assay_pair: str = "ptau217_vs_pet",
                          threshold: float | None = None,
                          missing_policy: str = "exclude") -> pd.Series:
    if assay_pair == "ptau217_vs_pet":
        column, default, positive = "PTAU217_bl", C.THRESHOLDS["p_tau217_pg_ml"], lambda x, t: x > t
    elif assay_pair == "abeta_ratio_vs_pet":
        column, default, positive = "AB_RATIO_bl", C.THRESHOLDS["abeta_ratio"], lambda x, t: x < t
    else:
        raise ValueError(f"未知 assay_pair: {assay_pair}")
    if missing_policy not in {"exclude", "missing_as_negative"}:
        raise ValueError(f"未知 missing_policy: {missing_policy}")
    values = pd.to_numeric(df[column], errors="coerce")
    status = positive(values, default if threshold is None else float(threshold)).astype("Int64")
    status[values.isna()] = pd.NA
    return status.fillna(0) if missing_policy == "missing_as_negative" else status


def add_canonical_group(df: pd.DataFrame, assay_pair: str = "ptau217_vs_pet",
                        threshold: float | None = None,
                        missing_policy: str = "exclude") -> pd.DataFrame:
    result = df.copy()
    result["PLASMA_STATUS"] = compute_plasma_status(result, assay_pair, threshold, missing_policy)
    result = result.dropna(subset=["PET_STATUS", "PLASMA_STATUS"]).copy()
    result["GROUP"] = result["PET_STATUS"].astype(int) * 2 + result["PLASMA_STATUS"].astype(int)
    return result


def _dated_table(key: str, value_cols: list[str]) -> pd.DataFrame:
    table = _read_csv(key).copy()
    date_col = C.FIELDS["date"][key]
    table[date_col] = pd.to_datetime(table[date_col], errors="coerce")
    table = table.dropna(subset=["RID", date_col])
    for col in value_cols:
        if col in table.columns and col not in {C.FIELDS.get("tracer"), "TRACER"}:
            table[col] = pd.to_numeric(table[col], errors="coerce")
    return table


def _nearest_measurement(index: pd.DataFrame, measures: pd.DataFrame, measure_date: str,
                         value_cols: list[str], max_days: int, prefix: str) -> pd.DataFrame:
    left = index[["RID", "INDEX_DATE"]].dropna().sort_values(["INDEX_DATE", "RID"])
    right = measures[["RID", measure_date] + value_cols].dropna(subset=[measure_date]).copy()
    right = right.rename(columns={measure_date: "MEASURE_DATE"}).sort_values(["MEASURE_DATE", "RID"])
    merged = pd.merge_asof(left, right, by="RID", left_on="INDEX_DATE", right_on="MEASURE_DATE",
                           direction="nearest", tolerance=pd.Timedelta(days=max_days))
    merged["TIME_GAP_DAYS"] = (merged["MEASURE_DATE"] - merged["INDEX_DATE"]).dt.days.abs()
    rename = {col: f"{prefix}_{col}" for col in value_cols}
    rename.update({"MEASURE_DATE": f"{prefix}_DATE", "TIME_GAP_DAYS": f"{prefix}_TIME_GAP_DAYS"})
    return merged.rename(columns=rename)


def build_aligned_cohort(window: str = "180d", assay_pair: str = "ptau217_vs_pet",
                         missing_policy: str = "exclude") -> pd.DataFrame:
    if window not in ALIGNMENT_WINDOWS:
        raise ValueError(f"window 必须是 {sorted(ALIGNMENT_WINDOWS)}")
    base = build_baseline()
    pet = _dated_table("pet", [C.FIELDS["pet_status"], C.FIELDS["tracer"], "SITEID"])
    pet[C.FIELDS["pet_status"]] = pd.to_numeric(pet[C.FIELDS["pet_status"]], errors="coerce")
    pet = pet.dropna(subset=[C.FIELDS["pet_status"]])
    pet = earliest_row(pet, C.FIELDS["date"]["pet"], [C.FIELDS["pet_status"]])
    pet = pet[["RID", C.FIELDS["date"]["pet"], "SITEID"]].rename(columns={
        C.FIELDS["date"]["pet"]: "INDEX_DATE", "SITEID": "PET_SITEID"})
    index = base.merge(pet, on="RID", how="inner")

    plasma_cols = [C.FIELDS["p_tau217"], C.FIELDS["abeta42_ab40"], "NfL_Q", "GFAP_Q"]
    aligned = _nearest_measurement(index, _dated_table("plasma", plasma_cols),
                                    C.FIELDS["date"]["plasma"], plasma_cols,
                                    ALIGNMENT_WINDOWS[window], "PLASMA")
    result = index.merge(aligned.drop(columns=["INDEX_DATE"]), on="RID", how="left")
    result["PTAU217_bl"] = clean_numeric(result[f"PLASMA_{C.FIELDS['p_tau217']}"])
    result["AB_RATIO_bl"] = clean_numeric(result[f"PLASMA_{C.FIELDS['abeta42_ab40']}"])
    result["GFAP_ALIGNED"] = clean_numeric(result["PLASMA_GFAP_Q"])
    result["NFL_ALIGNED"] = clean_numeric(result["PLASMA_NfL_Q"])
    result = add_canonical_group(result, assay_pair, missing_policy=missing_policy)

    tau_cols = [C.FIELDS["tau_suvr"], C.FIELDS["tracer"], C.FIELDS["qc"]]
    tau = _dated_table("tau", tau_cols)
    tau[C.FIELDS["qc"]] = pd.to_numeric(tau[C.FIELDS["qc"]], errors="coerce")
    tau = tau[(tau[C.FIELDS["tracer"]].astype(str).str.upper() == "FTP") & tau[C.FIELDS["qc"]].isin([1, 2])]
    tau_match = _nearest_measurement(index, tau, C.FIELDS["date"]["tau"], [C.FIELDS["tau_suvr"]], 365, "TAU")
    result = result.merge(tau_match.drop(columns=["INDEX_DATE"]), on="RID", how="left")
    result["TAU_SUVR_ALIGNED"] = pd.to_numeric(result[f"TAU_{C.FIELDS['tau_suvr']}"], errors="coerce")

    cdr = _dated_table("cdr", [C.FIELDS["cdrsb"]])
    adas = _dated_table("adas", [C.FIELDS["adas13"]])
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")
    result = result.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    result["D_CDRSB_yr"] = result["D_CDRSB"] / result["YRS_CDR"]
    result["D_ADAS13_yr"] = result["D_ADAS13"] / result["YRS_ADAS"]
    result["ALIGNMENT_WINDOW"] = window
    return result


def summarize_missingness(df: pd.DataFrame, columns: list[str]) -> dict:
    return {col: {"n": int(df[col].notna().sum()), "missing_rate": float(df[col].isna().mean())}
            for col in columns if col in df.columns}



def validate_raw_inputs() -> dict:
    required = {
        "ptdemog": ["RID", C.FIELDS["date"]["ptdemog"], C.FIELDS["pt_dob_yy"]],
        "cdr": ["RID", C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"]],
        "adas": ["RID", C.FIELDS["date"]["adas"], C.FIELDS["adas13"]],
        "plasma": ["RID", C.FIELDS["date"]["plasma"], C.FIELDS["p_tau217"], C.FIELDS["abeta42_ab40"]],
        "pet": ["RID", C.FIELDS["date"]["pet"], C.FIELDS["pet_status"]],
        "tau": ["RID", C.FIELDS["date"]["tau"], C.FIELDS["tau_suvr"], C.FIELDS["tracer"], C.FIELDS["qc"]],
    }
    report = {"valid": True, "tables": {}}
    for key, columns in required.items():
        path = C.RAW_DIR / C.FILES[key]
        item = {"path": str(path), "exists": path.exists(), "columns": {}, "n_rows": None}
        if not path.exists():
            report["valid"] = False
            report["tables"][key] = item
            continue
        table = pd.read_csv(path, nrows=5, low_memory=False)
        item["n_rows"] = int(sum(1 for _ in open(path, "rb")) - 1)
        item["columns"] = {col: col in table.columns for col in columns}
        if not all(item["columns"].values()):
            report["valid"] = False
        report["tables"][key] = item
    return report
