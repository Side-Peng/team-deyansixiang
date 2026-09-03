# -*- coding: utf-8 -*-
"""tau_analysis.py — 第五轮：tau PET 负担检验（PET+/Plasma− 组机制区分）

问题：PET+/Plasma−（551 人）认知轨迹接近双阳、合并症最高——
其"PET 阳性"是真实的 tau 病理连续谱，还是合并症/生物学变异驱动的假阳性？
tau PET（meta-temporal SUVR）是两者分界的关键证据。

设计：
- 主分析：TRACER=='FTP'(flortaucipir, 占92%) × qc_flag∈{1,2}（合格码），
  每 RID 取最早扫描的 META_TEMPORAL_SUVR（跨示踪剂 SUVR 不可比，警告字段已注明）
- 统计：四组 SUVR 分布（KW + 两两 MW，重点 PET+/Plasma− vs 双阴/双阳）
- tau 阳性率：1.22 / 1.28 / 1.36 三档敏感性（文献常用 FTP meta-temporal cutoff，
  非 ADNI 官方判定——官方只有 amyloid 的 AMYLOID_STATUS，tau 无官方阈值）
- 交叉验证：tau SUVR vs 血浆 p-tau217 的 Spearman 相关（不一致组内）
- 时间差敏感性：tau 扫描 vs amyloid 扫描日期差 >3 年者剔除重跑

用法：python scripts/tau_analysis.py
输出：data/processed/tau_report.txt
"""
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

import config as C

GROUPS = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}
TAU_CUTS = [1.22, 1.28, 1.36]  # 文献常用 FTP meta-temporal SUVR 阈值（敏感性三档）


def main():
    # ---------- 四组人群 ----------
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    sub = sub.dropna(subset=["PET_STATUS", "PLASMA_STATUS", "AGE"]).copy()
    sub["GROUP"] = (sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int))

    # ---------- tau 表：FTP + 合格 QC，每 RID 最早扫描 ----------
    tau = pd.read_csv(C.RAW_DIR / C.FILES["tau"], low_memory=False)
    n_all, n_rid = len(tau), tau.RID.nunique()
    tau = tau[tau["TRACER"] == "FTP"].copy()
    tau = tau[tau["qc_flag"].isin([1, 2])].copy()  # 合格码（2 为默认，1 为部分通过）
    tau[C.FIELDS["date"]["tau"]] = pd.to_datetime(tau[C.FIELDS["date"]["tau"]], errors="coerce")
    tau = tau.dropna(subset=[C.FIELDS["date"]["tau"], C.FIELDS["tau_suvr"]])
    tau = tau.sort_values(["RID", C.FIELDS["date"]["tau"]]).groupby("RID", as_index=False).first()
    tau = tau[["RID", C.FIELDS["date"]["tau"], C.FIELDS["tau_suvr"]]].rename(
        columns={C.FIELDS["date"]["tau"]: "TAU_DATE", C.FIELDS["tau_suvr"]: "TAU_SUVR"})

    # ---------- amyloid 基线日期（时间差用） ----------
    amy = pd.read_csv(C.RAW_DIR / C.FILES["pet"], low_memory=False)
    amy[C.FIELDS["date"]["pet"]] = pd.to_datetime(amy[C.FIELDS["date"]["pet"]], errors="coerce")
    amy = amy.dropna(subset=[C.FIELDS["date"]["pet"]]).sort_values(["RID", C.FIELDS["date"]["pet"]])
    amy = amy.groupby("RID", as_index=False).first()[["RID", C.FIELDS["date"]["pet"]]].rename(
        columns={C.FIELDS["date"]["pet"]: "AMY_DATE"})

    m = sub.merge(tau, on="RID", how="left").merge(amy, on="RID", how="left")
    m["TAU_DT"] = (m["TAU_DATE"] - m["AMY_DATE"]).dt.days

    lines = [f"tau 表原始: {n_all} 行 / {n_rid} RID → 筛选 FTP+QC合格+最早扫描: {len(tau)} 人"]

    # ---------- ① 可评估覆盖（选择偏倚记录） ----------
    lines.append("\n===== ① tau 可评估覆盖（FTP + QC 合格）=====")
    for g in sorted(GROUPS):
        d = m[m.GROUP == g]
        ev = d.TAU_SUVR.notna()
        lines.append(f"  {GROUPS[g]:<14} n={len(d):>4}  tau 可评估={ev.sum():>4}  ({ev.mean():.1%})")

    # ---------- ② SUVR 分布 + 两两检验 ----------
    lines.append("\n===== ② meta-temporal SUVR（每 RID 最早 FTP 扫描，中位数 [IQR]）=====")
    d = m.dropna(subset=["TAU_SUVR"])
    H, p = kruskal(*[d.loc[d.GROUP == g, "TAU_SUVR"].values for g in sorted(GROUPS)])
    lines.append(f"四组 KW: H={H:.1f}, p={p:.4g}")
    for g in sorted(GROUPS):
        v = d.loc[d.GROUP == g, "TAU_SUVR"]
        q = v.quantile([.5, .25, .75])
        lines.append(f"  {GROUPS[g]:<14} n={len(v):>4}  SUVR={q[.5]:.3f} [{q[.25]:.3f}, {q[.75]:.3f}]")

    pairs = [(2, 0, "PET+/Plasma− vs PET−/Plasma−"),
             (2, 3, "PET+/Plasma− vs PET+/Plasma+"),
             (1, 0, "PET−/Plasma+ vs PET−/Plasma−"),
             (1, 3, "PET−/Plasma+ vs PET+/Plasma+")]
    lines.append("两两 Mann-Whitney:")
    for g1, g2, name in pairs:
        a = d.loc[d.GROUP == g1, "TAU_SUVR"].values
        b = d.loc[d.GROUP == g2, "TAU_SUVR"].values
        p2 = mannwhitneyu(a, b, alternative="two-sided").pvalue
        lines.append(f"  {name}: p={p2:.4g}")

    # ---------- ③ tau 阳性率（三档敏感性） ----------
    lines.append("\n===== ③ tau 阳性率（FTP meta-temporal 文献常用阈值，三档敏感性；非 ADNI 官方判定）=====")
    hdr = f"  {'组':<14}" + "".join([f"{c:.2f}:{m2:>6}" for c, m2 in zip(TAU_CUTS, ["n", "n", "n"])])
    lines.append(f"  {'组':<14} {'n':>5}  " + "  ".join([f"≥{c:.2f} %" for c in TAU_CUTS]))
    for g in sorted(GROUPS):
        v = d.loc[d.GROUP == g, "TAU_SUVR"]
        pos = [f"{(v >= c).mean():.1%}" for c in TAU_CUTS]
        lines.append(f"  {GROUPS[g]:<14} {len(v):>5}  " + "  ".join(f"{p:>8}" for p in pos))

    # ---------- ④ 血浆 p-tau217 与 tau SUVR 一致性（不一致组内） ----------
    lines.append("\n===== ④ tau SUVR vs 血浆 p-tau217（Spearman，不一致组内）=====")
    for g in [1, 2]:
        dg = m[(m.GROUP == g)].dropna(subset=["TAU_SUVR", "PTAU217_bl"])
        if len(dg) >= 10:
            rho, ps = spearmanr(dg["TAU_SUVR"], dg["PTAU217_bl"])
            lines.append(f"  {GROUPS[g]:<14} n={len(dg):>3}  rho={rho:+.3f}, p={ps:.4g}")
        else:
            lines.append(f"  {GROUPS[g]}: 样本不足（{len(dg)}）")

    # ---------- ⑤ 时间差敏感性 ----------
    lines.append("\n===== ⑤ 时间差敏感性（tau vs amyloid 扫描）=====")
    lines.append(f"  时间差(天) 中位={m.TAU_DT.median():+.0f}，|差|>3年(1095天) 占比 {(m.TAU_DT.abs()>1095).mean():.1%}")
    d2 = d[d.TAU_DT.abs() <= 1095]
    H2, p2 = kruskal(*[d2.loc[d2.GROUP == g, "TAU_SUVR"].values for g in sorted(GROUPS) if len(d2[d2.GROUP == g]) >= 10])
    lines.append(f"  剔除 >3 年后（n={len(d2)}）四组 KW: H={H2:.1f}, p={p2:.4g}")
    for g in [0, 2, 3]:
        v = d2.loc[d2.GROUP == g, "TAU_SUVR"]
        if len(v):
            q = v.quantile([.5, .25, .75])
            lines.append(f"    {GROUPS[g]:<14} n={len(v):>4}  SUVR={q[.5]:.3f} [{q[.25]:.3f}, {q[.75]:.3f}]")

    # ---------- 结论核验 ----------
    lines.append("\n===== 机制判读（描述性，样本量注明）=====")
    med = {g: d.loc[d.GROUP == g, "TAU_SUVR"].median() for g in [0, 1, 2, 3]}
    lines.append(f"  双阴={med[0]:.2f} < P−/P+={med[1]:.2f} < PET+/P−={med[2]:.2f} ≈ 双阳={med[3]:.2f} ?")
    if med[2] > med[1] and med[2] <= med[3] * 1.05:
        lines.append("  → PET+/Plasma− tau 负担居中偏双阳：PET 阳性伴 tau 病理（真实连续谱）")
    elif med[2] <= med[0] * 1.1:
        lines.append("  → PET+/Plasma− tau 负担接近双阴：PET 阳性更可能是假阳性/非 tau 负荷")
    else:
        lines.append("  → 需结合 p 值与阈值阳性率综合判读")
    g2 = m[m.GROUP == 2]
    lines.append(f"  注意：PET+/Plasma− 仅 {g2.TAU_SUVR.notna().sum()}/{len(g2)} 可评估"
                 f"（{g2.TAU_SUVR.notna().mean():.1%}），选择偏倚是主要局限。")

    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "tau_report.txt").write_text(txt, encoding="utf-8")
    print("\n已保存: data/processed/tau_report.txt")


if __name__ == "__main__":
    main()
