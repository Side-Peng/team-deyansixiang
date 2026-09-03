# -*- coding: utf-8 -*-
"""tau_analysis_corrected.py — 修正口径下的 tau PET 机制区分（复赛 · bug 复核）

背景：初赛 tau_analysis.py 复用 subjects_wide.csv 的 PLASMA_STATUS 列，而该列由
preprocess.py 以 `PTAU217 > 0.183` 生成，pandas 中 `NaN > 阈值 = False`，导致
「未测血浆」被编码为「血浆阴性」，使 PET+/Plasma− 组从真 120 人膨胀到 551 人。

本脚本不改动初赛 tau_analysis.py，而是用正确口径（缺失血浆排除）重建四组后重跑
tau 分析，并同时跑「初赛口径（缺失=阴性）」做同口径对照，产出：
    data/processed/tau_report_corrected.txt   （修正口径完整报告 + 双口径对照表）

结论判读注意：tau 阳性阈值（1.22/1.28/1.36）为文献常用 FTP meta-temporal cutoff，
非 ADNI 官方判定；PET+/Plasma− 组 tau 可评估比例低（选择偏倚），须注明。
"""
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

import config as C
from environment import compute_plasma_status, GROUP_LABELS

TAU_CUTS = [1.22, 1.28, 1.36]


def load_tau_table():
    tau = pd.read_csv(C.RAW_DIR / C.FILES["tau"], low_memory=False)
    n_all, n_rid = len(tau), tau.RID.nunique()
    tau = tau[tau["TRACER"] == "FTP"].copy()
    tau = tau[tau["qc_flag"].isin([1, 2])].copy()
    tau[C.FIELDS["date"]["tau"]] = pd.to_datetime(tau[C.FIELDS["date"]["tau"]], errors="coerce")
    tau = tau.dropna(subset=[C.FIELDS["date"]["tau"], C.FIELDS["tau_suvr"]])
    tau = tau.sort_values(["RID", C.FIELDS["date"]["tau"]]).groupby("RID", as_index=False).first()
    tau = tau[["RID", C.FIELDS["date"]["tau"], C.FIELDS["tau_suvr"]]].rename(
        columns={C.FIELDS["date"]["tau"]: "TAU_DATE", C.FIELDS["tau_suvr"]: "TAU_SUVR"})
    return tau, n_all, n_rid


def load_amy_date():
    amy = pd.read_csv(C.RAW_DIR / C.FILES["pet"], low_memory=False)
    amy[C.FIELDS["date"]["pet"]] = pd.to_datetime(amy[C.FIELDS["date"]["pet"]], errors="coerce")
    amy = amy.dropna(subset=[C.FIELDS["date"]["pet"]]).sort_values(["RID", C.FIELDS["date"]["pet"]])
    amy = amy.groupby("RID", as_index=False).first()[["RID", C.FIELDS["date"]["pet"]]].rename(
        columns={C.FIELDS["date"]["pet"]: "AMY_DATE"})
    return amy


def build_groups(missing_policy: str) -> pd.DataFrame:
    """按给定 missing_policy 重建四组（含 AGE），返回带 GROUP 的基线帧。"""
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    plasma = compute_plasma_status(sub, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"], missing_policy)
    sub = sub.dropna(subset=["PET_STATUS", "AGE"]).copy()
    sub["PLASMA_STATUS"] = plasma.loc[sub.index]
    sub = sub.dropna(subset=["PET_STATUS", "PLASMA_STATUS"]).copy()
    sub["GROUP"] = (sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int))
    return sub


def analyze(sub: pd.DataFrame, tau: pd.DataFrame, amy: pd.DataFrame) -> dict:
    m = sub.merge(tau, on="RID", how="left").merge(amy, on="RID", how="left")
    m["TAU_DT"] = (m["TAU_DATE"] - m["AMY_DATE"]).dt.days
    n_by_group = {g: int((m.GROUP == g).sum()) for g in sorted(GROUP_LABELS)}

    # ① tau 可评估覆盖
    coverage = {}
    for g in sorted(GROUP_LABELS):
        d = m[m.GROUP == g]
        coverage[GROUP_LABELS[g]] = {"n": int(len(d)), "tau_evaluable": int(d.TAU_SUVR.notna().sum())}

    # ② SUVR 分布 + KW + 两两
    d = m.dropna(subset=["TAU_SUVR"])
    groups = [d.loc[d.GROUP == g, "TAU_SUVR"].values for g in sorted(GROUP_LABELS)]
    H, p = kruskal(*groups)
    med = {GROUP_LABELS[g]: float(d.loc[d.GROUP == g, "TAU_SUVR"].median()) for g in sorted(GROUP_LABELS)}
    n_eval = {GROUP_LABELS[g]: int((d.GROUP == g).sum()) for g in sorted(GROUP_LABELS)}
    pairs = [(2, 0, "PET+/Plasma− vs PET−/Plasma−"), (2, 3, "PET+/Plasma− vs PET+/Plasma+"),
             (1, 0, "PET−/Plasma+ vs PET−/Plasma−"), (1, 3, "PET−/Plasma+ vs PET+/Plasma+")]
    pairwise = {}
    for g1, g2, name in pairs:
        a = d.loc[d.GROUP == g1, "TAU_SUVR"].values
        b = d.loc[d.GROUP == g2, "TAU_SUVR"].values
        if len(a) and len(b):
            pairwise[name] = {"p": float(mannwhitneyu(a, b).pvalue), "n1": int(len(a)), "n2": int(len(b))}
        else:
            pairwise[name] = {"p": None, "n1": int(len(a)), "n2": int(len(b))}

    # ③ tau 阳性率（三档敏感性）
    pos_rates = {}
    for g in sorted(GROUP_LABELS):
        v = d.loc[d.GROUP == g, "TAU_SUVR"]
        pos_rates[GROUP_LABELS[g]] = {f"ge{c:.2f}": float((v >= c).mean()) for c in TAU_CUTS}

    # ④ 血浆 p-tau217 与 tau SUVR 一致性（不一致组内）
    spearman = {}
    for g in [1, 2]:
        dg = m[(m.GROUP == g)].dropna(subset=["TAU_SUVR", "PTAU217_bl"])
        if len(dg) >= 10:
            rho, ps = spearmanr(dg["TAU_SUVR"], dg["PTAU217_bl"])
            spearman[GROUP_LABELS[g]] = {"n": int(len(dg)), "rho": float(rho), "p": float(ps)}
        else:
            spearman[GROUP_LABELS[g]] = {"n": int(len(dg)), "rho": None, "p": None}

    # ⑤ 时间差敏感性
    d2 = d[d.TAU_DT.abs() <= 1095]
    gs2 = [d2.loc[d2.GROUP == g, "TAU_SUVR"].values for g in sorted(GROUP_LABELS) if len(d2[d2.GROUP == g]) >= 10]
    H2, p2 = (kruskal(*gs2) if len(gs2) >= 2 else (None, None))
    sens = {"n_after_exclude": int(len(d2)), "H": (float(H2) if H2 is not None else None),
            "p": (float(p2) if p2 is not None else None),
            "median_gt3yr": {GROUP_LABELS[g]: float(d2.loc[d2.GROUP == g, "TAU_SUVR"].median())
                             for g in [0, 2, 3] if len(d2[d2.GROUP == g])}}

    return {"n_by_group": n_by_group, "coverage": coverage,
            "KW": {"H": float(H), "p": float(p)}, "median_suvr": med, "n_evaluable": n_eval,
            "pairwise": pairwise, "pos_rates": pos_rates, "spearman": spearman, "sensitivity": sens}


def fmt_lines(label: str, res: dict) -> list:
    L = [f"\n{'='*64}", f"【{label}】", f"{'='*64}"]
    L.append(f"四组人数：{res['n_by_group']}")
    L.append("\n① tau 可评估覆盖：")
    for g, v in res["coverage"].items():
        L.append(f"  {g:<14} n={v['n']:>4}  tau可评估={v['tau_evaluable']:>4}")
    L.append(f"\n② meta-temporal SUVR（每 RID 最早 FTP 扫描，中位数）：")
    L.append(f"  KW H={res['KW']['H']:.1f}, p={res['KW']['p']:.4g}")
    for g in GROUP_LABELS.values():
        L.append(f"  {g:<14} n={res['n_evaluable'][g]:>4}  SUVR={res['median_suvr'][g]:.3f}")
    L.append("  两两 Mann-Whitney：")
    for name, v in res["pairwise"].items():
        L.append(f"    {name}: p={v['p'] if v['p'] is None else round(v['p'],4)} (n={v['n1']}/{v['n2']})")
    L.append("\n③ tau 阳性率（三档阈值，非 ADNI 官方判定）：")
    for g in GROUP_LABELS.values():
        r = res["pos_rates"][g]
        L.append(f"  {g:<14} " + "  ".join(f"≥{c:.2f}={r[f'ge{c:.2f}']:.1%}" for c in TAU_CUTS))
    L.append("\n④ tau SUVR vs 血浆 p-tau217（Spearman，不一致组内）：")
    for g in ["PET−/Plasma+", "PET+/Plasma−"]:
        s = res["spearman"][g]
        rho = "—" if s["rho"] is None else f"{s['rho']:+.3f}"
        pv = "—" if s["p"] is None else f"{s['p']:.4g}"
        L.append(f"  {g:<14} n={s['n']:>3}  rho={rho}  p={pv}")
    L.append("\n⑤ 时间差敏感性（剔除 |tau−amyloid|>3 年）：")
    L.append(f"  n={res['sensitivity']['n_after_exclude']}  KW H={res['sensitivity']['H']} p={res['sensitivity']['p']}")
    return L


def main():
    tau, n_all, n_rid = load_tau_table()
    amy = load_amy_date()

    corrected = analyze(build_groups("exclude"), tau, amy)
    legacy = analyze(build_groups("missing_as_negative"), tau, amy)

    lines = ["tau 表原始行数 / RID：%d / %d → FTP+QC合格+最早扫描：%d 人" % (n_all, n_rid, len(tau)),
             "\n# 修正口径（缺失血浆排除）"]
    lines += fmt_lines("修正口径 exclude", corrected)
    lines += fmt_lines("初赛口径 missing_as_negative（对照）", legacy)

    # 对照表：核心机制判读
    lines.append("\n" + "="*64)
    lines.append("核心对照：PET+/Plasma−（影像先行型）的 tau 负担在两个口径下是否一致")
    lines.append("="*64)
    for label, res in [("修正口径", corrected), ("初赛口径", legacy)]:
        m = res["median_suvr"]
        lines.append(f"\n{label}: 双阴={m['PET−/Plasma−']:.3f}  P−/P+={m['PET−/Plasma+']:.3f}  "
                     f"PET+/P−={m['PET+/Plasma−']:.3f}  双阳={m['PET+/Plasma+']:.3f}")
        pr = res["pairwise"]["PET+/Plasma− vs PET−/Plasma−"]
        p23 = res["pairwise"]["PET+/Plasma− vs PET+/Plasma+"]
        lines.append(f"  PET+/Plasma− vs 双阴 p={pr['p'] if pr['p'] is None else round(pr['p'],4)} (n={pr['n1']}/{pr['n2']})；"
                     f"vs 双阳 p={p23['p'] if p23['p'] is None else round(p23['p'],4)}")
    lines.append("\n判读提醒：PET+/Plasma− 组 tau 可评估比例低（选择偏倚）；tau 阈值为文献 cutoff 非 ADNI 官方。")

    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "tau_report_corrected.txt").write_text(txt, encoding="utf-8")
    print("\n已保存: data/processed/tau_report_corrected.txt")


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    main()
