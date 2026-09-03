# -*- coding: utf-8 -*-
"""preprocess.py — raw CSV → 基线（日期最早）合并 → 状态列 → subjects_wide.csv

用法：python scripts/preprocess.py
输出：data/processed/subjects_wide.csv, missing_report.csv
基线规则：每张表按 RID 取日期最早的非缺失行（各表 VISCODE 体系不统一，不用 'bl'）。
"""
import sys
import pandas as pd

import config as C

CLEAN_RULES = {  # 列 → 需排除的编码（ADNI 缺失编码：负值等）
    "pT217_F": lambda v: v <= 0,
    "AB42_AB40_F": lambda v: v <= 0,
}


def load(fname: str) -> pd.DataFrame:
    p = C.RAW_DIR / fname
    if not p.exists():
        sys.exit(f"[preprocess] 缺少文件: {p}")
    return pd.read_csv(p, low_memory=False)


def earliest(df: pd.DataFrame, date_col: str, cols: list) -> pd.DataFrame:
    """按 RID 取日期最早行，返回 [RID, *cols]。"""
    d = df.copy()
    d[date_col] = pd.to_datetime(d[date_col], errors="coerce")
    d = d.sort_values([C.FIELDS["rid"], date_col])
    # 取同一条最早源记录，避免 groupby().first() 跨访视逐列拼接。
    out = d.drop_duplicates(C.FIELDS["rid"], keep="first")
    keep = [C.FIELDS["rid"]] + [c for c in cols if c in d.columns]
    return out[keep]


def main():
    C.PROC_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 人口学（PTDEMOG 纵向，取每 RID 最早行）----
    demo = load(C.FILES["ptdemog"])
    d0 = earliest(demo, C.FIELDS["date"]["ptdemog"],
                  [C.FIELDS["pt_dob_yy"], C.FIELDS["gender"], C.FIELDS["educat"], "VISDATE"])
    # PTDOBYY 是出生年份 float（如 1931.0）：不能喂给 pd.to_datetime（会被当作纳秒时间戳 → 1970 年），
    # 年份数值直接相减即可。
    dob_yy = pd.to_numeric(d0[C.FIELDS["pt_dob_yy"]], errors="coerce")
    d0["AGE"] = pd.to_datetime(d0["VISDATE"]).dt.year - dob_yy
    out = d0[[C.FIELDS["rid"], "AGE", C.FIELDS["gender"], C.FIELDS["educat"]]].rename(
        columns={C.FIELDS["rid"]: "RID", C.FIELDS["gender"]: "GENDER",
                 C.FIELDS["educat"]: "EDUCAT"})
    out["GENDER"] = C.normalize_gender(out["GENDER"])  # ADNIMERGE2 字符串 → 数字编码

    # ---- 认知与诊断（各表取日期最早行）----
    for tbl, val_col, out_col in [
        ("cdr", C.FIELDS["cdrsb"], "CDRSB_bl"),
        ("mmse", C.FIELDS["mmse"], "MMSE_bl"),
        ("adas", C.FIELDS["adas13"], "ADAS13_bl"),
        ("dxsum", C.FIELDS["dx"], "DX_bl"),
    ]:
        t = load(C.FILES[tbl])
        e = earliest(t, C.FIELDS["date"][tbl], [val_col])
        if tbl == "dxsum":
            e[out_col] = C.normalize_dx(e[val_col])  # ADNIMERGE2 'CN'/'MCI'/'Dementia' → 1/2/3
        else:
            e[out_col] = pd.to_numeric(e[val_col], errors="coerce")
        out = out.merge(e[[C.FIELDS["rid"], out_col]], on=C.FIELDS["rid"], how="left")

    # ---- 血浆（每 RID 日期最早行；清洗负值缺失编码）----
    pl = load(C.FILES["plasma"])
    pl_cols = [C.FIELDS["p_tau217"], C.FIELDS["abeta42_ab40"]]
    for c in pl_cols:
        pl.loc[CLEAN_RULES[c](pl[c]), c] = pd.NA
    e = earliest(pl, C.FIELDS["date"]["plasma"], pl_cols)
    out = out.merge(e[[C.FIELDS["rid"]] + pl_cols].rename(
        columns={C.FIELDS["p_tau217"]: "PTAU217_bl", C.FIELDS["abeta42_ab40"]: "AB_RATIO_bl"}),
        on=C.FIELDS["rid"], how="left")

    # ---- Aβ PET（每 RID 日期最早、状态非缺失）----
    pet = load(C.FILES["pet"])
    pet = pet[pet[C.FIELDS["pet_status"]].notna()].copy()
    e = earliest(pet, C.FIELDS["date"]["pet"], [C.FIELDS["pet_status"], C.FIELDS["centiloids"], C.FIELDS["tracer"]])
    out = out.merge(e[[C.FIELDS["rid"], C.FIELDS["pet_status"], C.FIELDS["centiloids"], C.FIELDS["tracer"]]].rename(
        columns={C.FIELDS["pet_status"]: "PET_STATUS", C.FIELDS["centiloids"]: "CENTILOIDS_bl",
                 C.FIELDS["tracer"]: "PET_TRACER"}),
        on=C.FIELDS["rid"], how="left")
    out["PET_STATUS"] = pd.to_numeric(out["PET_STATUS"], errors="coerce")

    # ---- 状态卡片：阈值未配置则不生成血浆状态（不编造阈值）----
    if C.THRESHOLDS["p_tau217_pg_ml"]:
        out["PLASMA_STATUS"] = (out["PTAU217_bl"] > C.THRESHOLDS["p_tau217_pg_ml"]).astype("Int64")
        out.loc[out["PTAU217_bl"].isna(), "PLASMA_STATUS"] = pd.NA
    else:
        print("[preprocess] ⚠️ 未生成 PLASMA_STATUS：先运行 calibrate.py 并将 Youden 阈值填入 config.THRESHOLDS")

    out.to_csv(C.PROC_DIR / "subjects_wide.csv", index=False)
    miss = out.isna().mean().round(3).to_frame("missing_rate")
    miss.to_csv(C.PROC_DIR / "missing_report.csv")
    print(f"[preprocess] 完成：基线 {len(out)} 人，关键字段缺失率：")
    print(miss.loc[["AGE", "CDRSB_bl", "MMSE_bl", "ADAS13_bl", "PTAU217_bl", "AB_RATIO_bl", "PET_STATUS"]].to_string())


if __name__ == "__main__":
    main()
