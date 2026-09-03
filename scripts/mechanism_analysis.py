# -*- coding: utf-8 -*-
"""Mechanism profiling for biomarker discordance."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.linear_model import HuberRegressor

import config as C
from data_layer import GROUP_LABELS, build_aligned_cohort, summarize_missingness

MECHANISM_GROUPS = {"PET−/Plasma+": 1, "PET+/Plasma−": 2}
MARKER_COLUMNS = {"tau_pet": "TAU_SUVR_ALIGNED", "GFAP": "GFAP_ALIGNED", "NfL": "NFL_ALIGNED"}
OUTCOME_COLUMNS = {"ADAS13": "D_ADAS13_yr", "CDRSB": "D_CDRSB_yr"}


def _finite(values) -> np.ndarray:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    return arr[np.isfinite(arr)]


def _safe_median(values):
    arr = _finite(values)
    return float(np.median(arr)) if len(arr) else None


def _effect(target, reference, seed=2026, n_boot=400):
    a, b = _finite(target), _finite(reference)
    result = {"n_target": int(len(a)), "n_reference": int(len(b)), "target_median": _safe_median(a),
              "reference_median": _safe_median(b), "median_difference": None,
              "mannwhitney_p": None, "cohens_d": None, "ci95": [None, None]}
    if not len(a) or not len(b):
        return result
    result["median_difference"] = float(np.median(a) - np.median(b))
    result["mannwhitney_p"] = float(sps.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    pooled = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) /
                     max(len(a) + len(b) - 2, 1))
    result["cohens_d"] = float((np.mean(a) - np.mean(b)) / pooled) if pooled > 0 else 0.0
    rng = np.random.default_rng(seed)
    diffs = [float(np.median(rng.choice(a, len(a), replace=True)) -
                   np.median(rng.choice(b, len(b), replace=True))) for _ in range(n_boot)]
    result["ci95"] = [float(x) for x in np.percentile(diffs, [2.5, 97.5])]
    return result


def _bh_adjust(p_values: dict) -> dict:
    valid = [(key, value) for key, value in p_values.items() if value is not None and np.isfinite(value)]
    if not valid:
        return {}
    valid.sort(key=lambda item: item[1])
    n = len(valid)
    adjusted = {}
    running = 1.0
    for rank, (key, value) in reversed(list(enumerate(valid, start=1))):
        running = min(running, value * n / rank)
        adjusted[key] = float(min(running, 1.0))
    return adjusted


def _regression_effect(frame: pd.DataFrame, outcome: str, target_group: int) -> dict:
    columns = [outcome, "GROUP", "AGE", "GENDER", "EDUCAT"]
    data = frame.dropna(subset=[col for col in columns if col in frame.columns]).copy()
    if len(data) < 20:
        return {"n": int(len(data)), "adjusted_difference": None, "r2": None}
    x = pd.DataFrame({"target": (data["GROUP"] == target_group).astype(int),
                      "AGE": pd.to_numeric(data["AGE"], errors="coerce"),
                      "GENDER": data["GENDER"].map({1: 1, 2: 0}),
                      "EDUCAT": pd.to_numeric(data["EDUCAT"], errors="coerce")}, index=data.index).dropna()
    if len(x) < 20:
        return {"n": int(len(x)), "adjusted_difference": None, "r2": None}
    y = pd.to_numeric(data.loc[x.index, outcome], errors="coerce")
    model = HuberRegressor(epsilon=1.35, alpha=0.0, max_iter=1000).fit(x, y)
    return {"n": int(len(x)), "model": "HuberRegressor", "adjusted_difference": float(model.coef_[0]),
            "r2": float(model.score(x, y))}


def _spearman_association(frame: pd.DataFrame, left: str, right: str) -> dict:
    data = frame[[left, right]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < 10:
        return {"n": int(len(data)), "rho": None, "p": None}
    rho, p = sps.spearmanr(data[left], data[right])
    return {"n": int(len(data)), "rho": float(rho), "p": float(p)}


def _mechanism_ranking(evidence: dict) -> list[dict]:
    lag = evidence["lag_evidence"]
    injury = evidence["injury_evidence"]
    tau = evidence["tau_evidence"]
    trajectory = evidence["trajectory_evidence"]
    scores = {"temporal_lag": 0.0, "non_ad_injury": 0.0, "measurement_noise": 0.0}
    reasons = {key: [] for key in scores}
    tau_diff = tau.get("effect", {}).get("median_difference")
    plasma_diff = lag.get("plasma_tau_vs_reference", {}).get("median_difference")
    injury_diffs = [injury.get(marker, {}).get("effect", {}).get("median_difference") for marker in ["GFAP", "NfL"]]
    if tau_diff is not None and tau_diff > 0 and (plasma_diff is None or plasma_diff <= 0):
        scores["temporal_lag"] += 1.0
        reasons["temporal_lag"].append("tau PET 升高而 plasma p-tau217 未同步升高")
    elif tau_diff is not None and tau_diff <= 0:
        scores["measurement_noise"] += 0.5
        reasons["measurement_noise"].append("tau PET 未显示同步升高")
    elif tau_diff is not None and plasma_diff is not None:
        reasons["temporal_lag"].append("tau PET 与 plasma p-tau217 同向，暂不支持明显滞后")
    if any(value is not None and value > 0 for value in injury_diffs):
        scores["non_ad_injury"] += 1.0
        reasons["non_ad_injury"].append("GFAP/NfL 至少一个相对参考组升高")
    if lag.get("plasma_tau_vs_reference", {}).get("mannwhitney_p") is not None:
        scores["temporal_lag"] += 0.25
    if trajectory.get("adjusted_difference") is not None and abs(trajectory["adjusted_difference"]) < 0.2:
        scores["measurement_noise"] += 0.5
        reasons["measurement_noise"].append("认知轨迹校正后差异较小")
    ranked = sorted(scores, key=scores.get, reverse=True)
    return [{"mechanism": key, "score": round(scores[key], 3), "reasons": reasons[key],
             "interpretation": "候选机制，需外部验证"} for key in ranked]


def profile_mechanism(frame: pd.DataFrame, target_group: str = "PET−/Plasma+",
                      markers: list[str] | None = None, outcome: str = "ADAS13",
                      alignment_window: str = "180d", seed: int = 2026) -> dict:
    if target_group not in MECHANISM_GROUPS:
        raise ValueError(f"target_group 必须是 {sorted(MECHANISM_GROUPS)}")
    markers = markers or ["tau_pet", "GFAP", "NfL"]
    unknown = sorted(set(markers) - set(MARKER_COLUMNS))
    if unknown:
        raise ValueError(f"未知 markers: {unknown}")
    if outcome not in OUTCOME_COLUMNS:
        raise ValueError(f"未知 outcome: {outcome}")
    target = MECHANISM_GROUPS[target_group]
    reference = frame[frame["GROUP"] == 0]
    target_frame = frame[frame["GROUP"] == target]
    effects = {}
    p_values = {}
    for marker in markers:
        column = MARKER_COLUMNS[marker]
        effect = _effect(target_frame[column], reference[column], seed=seed)
        effects[marker] = {"column": column, "effect": effect,
                           "missingness": {"target": float(target_frame[column].isna().mean()),
                                           "reference": float(reference[column].isna().mean())}}
        p_values[f"{marker}_vs_reference"] = effect["mannwhitney_p"]

    outcome_column = OUTCOME_COLUMNS[outcome]
    trajectory_effect = _effect(target_frame[outcome_column], reference[outcome_column], seed=seed + 1)
    trajectory_adjusted = _regression_effect(frame, outcome_column, target)
    p_values["trajectory"] = trajectory_effect["mannwhitney_p"]
    tau_column = MARKER_COLUMNS["tau_pet"]
    plasma_tau_column = "PTAU217_bl"
    plasma_tau_effect = _effect(target_frame[plasma_tau_column], reference[plasma_tau_column], seed=seed + 2)
    tau_lag_effect = _effect(target_frame[tau_column], reference[tau_column], seed=seed + 3)
    p_values["plasma_tau_vs_reference"] = plasma_tau_effect["mannwhitney_p"]
    p_values["tau_pet_vs_reference"] = tau_lag_effect["mannwhitney_p"]

    evidence = {
        "target_group": target_group,
        "target_group_code": target,
        "reference_group": GROUP_LABELS[0],
        "alignment_window": alignment_window,
        "sample": {"n_total": int(len(frame)), "n_target": int(len(target_frame)),
                    "n_reference": int(len(reference))},
        "lag_evidence": {"plasma_tau_vs_reference": plasma_tau_effect,
                         "tau_pet_vs_reference": tau_lag_effect,
                         "time_gap_days": {"plasma_median": _safe_median(target_frame.get("PLASMA_TIME_GAP_DAYS", [])),
                                            "tau_median": _safe_median(target_frame.get("TAU_TIME_GAP_DAYS", []))}},
        "injury_evidence": {marker: effects[marker] for marker in ["GFAP", "NfL"] if marker in effects},
        "tau_evidence": effects.get("tau_pet", {"effect": _effect([], [])}),
        "trajectory_evidence": {"outcome": outcome, "effect": trajectory_effect,
                                 "adjusted": trajectory_adjusted},
        "missingness_evidence": summarize_missingness(frame, [MARKER_COLUMNS[m] for m in markers] +
                                                       [outcome_column, "PTAU217_bl"]),
        "continuous_value_evidence": {
            "pTau217_vs_tau_pet": _spearman_association(target_frame, "PTAU217_bl", tau_column),
            "pTau217_vs_GFAP": _spearman_association(target_frame, "PTAU217_bl", "GFAP_ALIGNED"),
            "pTau217_vs_NfL": _spearman_association(target_frame, "PTAU217_bl", "NFL_ALIGNED"),
        },
    }
    evidence["multiple_testing"] = {"method": "Benjamini-Hochberg", "q_values": _bh_adjust(p_values)}
    evidence["candidate_mechanism_ranking"] = _mechanism_ranking(evidence)
    evidence["negative_evidence"] = [
        "机制评分仅用于排序，不等同于因果证明。",
        "时间对齐窗口和 tau PET 选择偏倚仍可能影响方向。",
    ]
    evidence["next_best_action"] = {
        "action": "profile_mechanism",
        "params": {"target_group": "PET+/Plasma−" if target_group == "PET−/Plasma+" else "PET−/Plasma+",
                   "markers": markers, "alignment_window": alignment_window, "outcome": outcome},
    }
    return evidence


def run_mechanism_analysis(window: str = "180d", output: Path | None = None) -> dict:
    frame = build_aligned_cohort(window=window, missing_policy="exclude")
    from data_layer import build_baseline, add_canonical_group
    canonical = add_canonical_group(build_baseline(), missing_policy="exclude")
    report = {"window": window, "n_aligned": int(len(frame)),
              "canonical_n": int(len(canonical)),
              "canonical_group_counts": {GROUP_LABELS[g]: int((canonical.GROUP == g).sum()) for g in GROUP_LABELS},
              "group_counts": {GROUP_LABELS[g]: int((frame.GROUP == g).sum()) for g in GROUP_LABELS},
              "profiles": {group: profile_mechanism(frame, group, alignment_window=window)
                           for group in MECHANISM_GROUPS}}
    if output is None:
        output = C.PROC_DIR / "mechanism_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    result = run_mechanism_analysis(sys.argv[1] if len(sys.argv) > 1 else "180d")
    print(json.dumps(result, ensure_ascii=False, indent=2)[:6000])
