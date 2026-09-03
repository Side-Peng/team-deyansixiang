# -*- coding: utf-8 -*-
"""Patient-level site-held-out validation for mechanism profiles."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import config as C
from data_layer import GROUP_LABELS, build_aligned_cohort
from mechanism_analysis import MECHANISM_GROUPS, profile_mechanism


def _site_split(frame):
    if "PET_SITEID" not in frame.columns or frame["PET_SITEID"].nunique(dropna=True) < 2:
        raise ValueError("site-held-out validation 至少需要两个非空 PET_SITEID")
    sites = sorted(frame["PET_SITEID"].dropna().astype(str).unique())
    holdout = {site for index, site in enumerate(sites) if index % 5 == 0}
    if len(holdout) == len(sites):
        holdout = set(sites[-1:])
    train = frame[~frame["PET_SITEID"].astype(str).isin(holdout)].copy()
    valid = frame[frame["PET_SITEID"].astype(str).isin(holdout)].copy()
    return train, valid, sorted(holdout)


def _ranking_names(profile):
    return [item["mechanism"] for item in profile.get("candidate_mechanism_ranking", [])]


def run_validation(window: str = "180d", output: Path | None = None) -> dict:
    frame = build_aligned_cohort(window=window, missing_policy="exclude")
    train, valid, holdout_sites = _site_split(frame)
    report = {"window": window, "split": {"unit": "RID", "strategy": "PET_SITEID holdout",
              "holdout_sites": holdout_sites, "n_train": int(len(train)), "n_validation": int(len(valid))},
              "fixed_thresholds": True, "group_counts": {}, "profiles": {}, "ranking_reproduction": {}}
    for group_name in MECHANISM_GROUPS:
        train_profile = profile_mechanism(train, group_name, alignment_window=window)
        valid_profile = profile_mechanism(valid, group_name, alignment_window=window)
        report["profiles"][group_name] = {"train": train_profile, "validation": valid_profile}
        train_counts = {GROUP_LABELS[g]: int((train.GROUP == g).sum()) for g in GROUP_LABELS}
        valid_counts = {GROUP_LABELS[g]: int((valid.GROUP == g).sum()) for g in GROUP_LABELS}
        report["group_counts"][group_name] = {"train": train_counts, "validation": valid_counts}
        train_rank = _ranking_names(train_profile)
        valid_rank = _ranking_names(valid_profile)
        overlap = len(set(train_rank[:2]) & set(valid_rank[:2])) / 2 if valid_rank else 0.0
        report["ranking_reproduction"][group_name] = {"train_top2": train_rank[:2],
                                                        "validation_top2": valid_rank[:2],
                                                        "top2_overlap": float(overlap)}
    report["limitations"] = [
        "site-held-out 是内部验证，不等同于独立外部队列。",
        "阈值固定，但机制排名仍属于假设生成。",
        "同一 RID 仅出现在一个 split；alignment 后以 index PET site 作为分层单位。",
    ]
    if output is None:
        output = C.PROC_DIR / "validation_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print(json.dumps(run_validation(sys.argv[1] if len(sys.argv) > 1 else "180d"), ensure_ascii=False, indent=2)[:6000])
