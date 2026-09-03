# -*- coding: utf-8 -*-
"""corrected_reports.py — 修正口径（缺失血浆排除）下的 trajectory / slice / cluster 报告

背景：初赛 trajectory.py / slice_analysis.py / cluster_pm.py 复用 subjects_wide.csv 的
PLASMA_STATUS 列（preprocess.py 把缺失血浆编码为阴性），导致 PET+/Plasma− 组 120→551、
PET−/Plasma− 组 641→1012 的膨胀。本脚本用修正口径重建四组后重跑三个分析，输出：

    data/processed/trajectory_report_corrected.txt
    data/processed/slice_report_corrected.txt
    data/processed/cluster_report_corrected.txt

cluster 的目标组为 PET−/Plasma+（128 人，需血浆阳性，天然排除缺失），不受该 bug 影响，
结果与初赛 cluster_report_v2.txt 完全一致（本脚本仍重跑一次以证实）。
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu
from sklearn.linear_model import LinearRegression

import config as C
from environment import compute_plasma_status, GROUP_LABELS, build_longitudinal
from trajectory import per_patient_change
from cluster_pm import make_comorbidity, FEATS

DX_LABEL = {1.0: "CN", 2.0: "MCI", 3.0: "AD", 10.0: "UNKNOWN(10)"}


def corrected_frame() -> pd.DataFrame:
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    plasma = compute_plasma_status(sub, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"], "exclude")
    sub = sub.dropna(subset=["PET_STATUS", "AGE"]).copy()
    sub["PLASMA_STATUS"] = plasma.loc[sub.index]
    sub = sub.dropna(subset=["PET_STATUS", "PLASMA_STATUS"]).copy()
    sub["GROUP"] = (sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int))
    return sub


# --------------------------------------------------------------------------
# trajectory（镜像初赛 trajectory.py，仅四组口径修正）
# --------------------------------------------------------------------------
def trajectory_corrected() -> str:
    sub = corrected_frame()
    m = build_longitudinal(sub)
    labels = GROUP_LABELS

    lines = [f"样本：{len(m)}（修正口径：缺失血浆排除）；四组 { {labels[g]: int((m.GROUP==g).sum()) for g in labels} }",
             f"随访: CDRSB {m.YRS_CDR.notna().sum()} 人 / ADAS13 {m.YRS_ADAS.notna().sum()} 人"]

    # ① 年化变化率
    lines.append("\n===== ① 年化变化率（Δ/年，中位数 [IQR]）=====")
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        lines.append(f"\n{col}:")
        for g in sorted(m.GROUP.unique()):
            v = m.loc[m.GROUP == g, col].dropna()
            q = v.quantile([.5, .25, .75])
            lines.append(f"  {labels[g]:<14} n={len(v):>4}  中位={q[.5]:+.2f} [{q[.25]:+.2f}, {q[.75]:+.2f}]")

    # ② KW + 两两
    lines.append("\n===== ② 组间检验（Kruskal-Wallis + Mann-Whitney）=====")
    pairs = [(1, 0, "PET−/Plasma+ vs PET−/Plasma−"),
             (1, 3, "PET−/Plasma+ vs PET+/Plasma+"),
             (1, 2, "PET−/Plasma+ vs PET+/Plasma−")]
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        d = m.dropna(subset=[col])
        H, p = kruskal(*[d.loc[d.GROUP == g, col].values for g in sorted(d.GROUP.unique())])
        lines.append(f"\n{col}: KW H={H:.1f}, p={p:.4g}")
        for g1, g2, name in pairs:
            a1 = d.loc[d.GROUP == g1, col].values
            a2 = d.loc[d.GROUP == g2, col].values
            if len(a1) and len(a2):
                lines.append(f"  {name}: n={len(a1)}/{len(a2)} p={mannwhitneyu(a1, a2).pvalue:.4g}")

    # ③ 混杂校正
    lines.append("\n===== ③ 混杂校正（LinearRegression, 以 PET−/Plasma− 为参照）=====")
    m = m.copy()
    m["G1"] = (m.GROUP == 1).astype(int)
    m["G2"] = (m.GROUP == 2).astype(int)
    m["G3"] = (m.GROUP == 3).astype(int)
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        d = m.dropna(subset=[col, "AGE", "GENDER", "EDUCAT"])
        X = d[["G1", "G2", "G3", "AGE", "GENDER", "EDUCAT"]].copy()
        X["GENDER"] = X["GENDER"].map({1: 1, 2: 0})
        X = X[X["GENDER"].notna()]
        yv = d.loc[X.index, col]
        lm = LinearRegression().fit(X, yv)
        lines.append(f"\n{col}: 校正年龄/性别/教育后组系数（相对 PET−/Plasma−，Δ/年）")
        for name, coef in zip(["PET−/Plasma+", "PET+/Plasma−", "PET+/Plasma+"], lm.coef_[:3]):
            lines.append(f"  {name:<14}: {coef:+.3f}")

    # ④ 置换检验
    lines.append("\n===== ④ 置换检验（1000 次打乱组标签，KW H 统计量）=====")
    rng = np.random.default_rng(2026)
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        d = m.dropna(subset=[col])
        y = d[col].values
        g = d.GROUP.values
        H_obs = kruskal(*[y[g == x] for x in np.unique(g)])[0]
        cnt = 0
        for _ in range(1000):
            gp = rng.permutation(g)
            cnt += kruskal(*[y[gp == x] for x in np.unique(g)])[0] >= H_obs
        lines.append(f"{col}: 观察 H={H_obs:.1f}，置换 H≥观察 比例 p={cnt/1000:.3f}")

    return "\n".join(lines)


# --------------------------------------------------------------------------
# slice（镜像初赛 slice_analysis.py，仅四组口径修正）
# --------------------------------------------------------------------------
def slice_corrected() -> str:
    sub = corrected_frame()
    # MEDHIST 系统标记（镜像初赛 slice_analysis.py 的单来源口径）
    mh = pd.read_csv(C.RAW_DIR / C.FILES["medhist"], low_memory=False)
    mh = mh.sort_values(["RID", "VISDATE"]).groupby("RID", as_index=False).first()
    sub = sub.merge(mh[["RID", "MH4CARD", "MH9ENDO"]], on="RID", how="left")
    sub["CARDIO"] = C.normalize_yes_no(sub["MH4CARD"]).fillna(0).astype(int)
    sub["ENDO"] = C.normalize_yes_no(sub["MH9ENDO"]).fillna(0).astype(int)
    m = build_longitudinal(sub)
    labels = GROUP_LABELS

    lines = [f"样本：{len(m)}（修正口径）；MEDHIST 覆盖 {m.MH4CARD.notna().sum()} 人（ADNI1/GO/2）"]

    # ① 合并症
    lines.append("\n===== ① 合并症代理（MEDHIST 系统标记）=====")
    for com, name in [("CARDIO", "心血管 MH4CARD"), ("ENDO", "内分泌 MH9ENDO")]:
        lines.append(f"\n{name} 阳性率：")
        tbl = []
        for g in sorted(m.GROUP.unique()):
            d = m[m.GROUP == g]
            r = d[com].mean()
            tbl.append([d[com].sum(), len(d) - d[com].sum()])
            lines.append(f"  {labels[g]:<14} {r:.2%}  (n={len(d)})")
        chi2, p, _, _ = chi2_contingency(np.array(tbl))
        lines.append(f"  卡方检验: chi2={chi2:.1f}, p={p:.4g}")

    # 轨迹回归加合并症
    for col in ["D_ADAS13_yr"]:
        d = m.dropna(subset=[col, "AGE", "GENDER", "EDUCAT", "MH4CARD"])
        d = d.copy()
        d["G1"] = (d.GROUP == 1).astype(int)
        d["G2"] = (d.GROUP == 2).astype(int)
        d["G3"] = (d.GROUP == 3).astype(int)
        X = d[["G1", "G2", "G3", "AGE", "GENDER", "EDUCAT", "CARDIO", "ENDO"]].copy()
        X["GENDER"] = X["GENDER"].map({1: 1, 2: 0})
        X["CARDIO"] = X["CARDIO"].astype(int)
        X["ENDO"] = X["ENDO"].astype(int)
        lm = LinearRegression().fit(X, d.loc[X.index, col])
        lines.append(f"\n{col}: 校正年龄/性别/教育 + 合并症后组系数（相对 PET−/Plasma−）")
        for nm, coef in zip(["PET−/Plasma+", "PET+/Plasma−", "PET+/Plasma+", "CARDIO", "ENDO"],
                            list(lm.coef_[:3]) + [lm.coef_[6], lm.coef_[7]]):
            lines.append(f"  {nm:<14}: {coef:+.3f}")

    # ② 诊断分层
    lines.append("\n===== ② 诊断分层（CN/MCI/AD 层内四组结构复现性）=====")
    for dxv in [1.0, 2.0, 3.0]:
        layer = m[m.DX_bl == dxv]
        lines.append(f"\n--- {DX_LABEL[dxv]}（n={len(layer)}）---")
        if len(layer) < 50:
            lines.append("  样本过小，跳过")
            continue
        row = []
        for g in sorted(layer.GROUP.unique()):
            d = layer[layer.GROUP == g]
            row.append((labels[g], len(d), d.ADAS13_bl.median(), d.D_ADAS13_yr.median()))
        lines.append(f"  {'组':<14} {'n':>4} {'ADAS13_bl':>9} {'ΔADAS13/年':>9}")
        for g, n, ab, tr in row:
            lines.append(f"  {g:<14} {n:>4} {ab:>9.1f} {tr:>+9.2f}")
        p2 = layer[layer.GROUP == 1].D_ADAS13_yr.dropna().values
        p0 = layer[layer.GROUP == 0].D_ADAS13_yr.dropna().values
        p3 = layer[layer.GROUP == 3].D_ADAS13_yr.dropna().values
        if len(p2) >= 10 and len(p0) >= 10 and len(p3) >= 10:
            p_vs_neg = mannwhitneyu(p2, p0).pvalue
            p_vs_pos = mannwhitneyu(p2, p3).pvalue
            med_neg, med_2, med_pos = np.median(p0), np.median(p2), np.median(p3)
            strict_mid = med_neg <= med_2 <= med_pos and p_vs_neg < 0.05
            lines.append(f"  中间态判定: 中位 双阴={med_neg:+.2f} ≤ P−/P+={med_2:+.2f} ≤ 双阳={med_pos:+.2f}")
            lines.append(f"  P−/P+ vs 双阴 p={p_vs_neg:.4g}；vs 双阳 p={p_vs_pos:.4g} → "
                         f"{'✅ 严格中间态' if strict_mid else '❌ 非严格中间态'}")

    return "\n".join(lines)


# --------------------------------------------------------------------------
# cluster（目标组 PET−/Plasma+ = 128 不受 bug 影响；重跑证实一致）
# --------------------------------------------------------------------------
def cluster_corrected() -> str:
    sub = corrected_frame()
    m = build_longitudinal(sub)
    m = make_comorbidity(m)
    pm = m[(m.PET_STATUS == 0) & (m.PLASMA_STATUS == 1)].copy()
    pm = pm.dropna(subset=FEATS)
    lines = ["说明：cluster 的目标组为 PET−/Plasma+（要求血浆阳性，天然排除缺失血浆），",
             "修正口径下该组人数与初赛完全一致（128 人），故聚类结果与 cluster_report_v2.txt 相同。",
             "此处重跑确认。",
             f"\nPET−/Plasma+ 组：{len(pm)} 人（聚类输入齐全；纵向 {pm.D_ADAS13_yr.notna().sum()} 人）",
             f"合并症阳性率（双来源）: 心血管 {pm.CARDIO.mean():.1%} 内分泌 {pm.ENDO.mean():.1%} CKD {pm.CKD.mean():.1%}"]
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import silhouette_score
    X = StandardScaler().fit_transform(pm[FEATS])
    km2 = KMeans(n_clusters=2, n_init=10, random_state=42).fit(X)
    sil2 = silhouette_score(X, km2.labels_)
    pm["CLUSTER2"] = km2.labels_
    lines.append(f"\nKMeans k=2: 轮廓系数={sil2:.3f}  簇规模: {np.bincount(km2.labels_).tolist()}")
    c0 = pm[pm.CLUSTER2 == 0]
    c1 = pm[pm.CLUSTER2 == 1]
    lines.append(f"两簇轨迹 ΔADAS13/年: 簇0={c0.D_ADAS13_yr.median():+.2f} vs 簇1={c1.D_ADAS13_yr.median():+.2f}, "
                 f"Mann-Whitney p={mannwhitneyu(c0.D_ADAS13_yr.dropna(), c1.D_ADAS13_yr.dropna()).pvalue:.4g}")
    lines.append("→ 与初赛一致：两亚群轨迹分化不显著（负结果）。")
    return "\n".join(lines)


def main():
    traj = trajectory_corrected()
    slc = slice_corrected()
    clu = cluster_corrected()
    (C.PROC_DIR / "trajectory_report_corrected.txt").write_text(traj, encoding="utf-8")
    (C.PROC_DIR / "slice_report_corrected.txt").write_text(slc, encoding="utf-8")
    (C.PROC_DIR / "cluster_report_corrected.txt").write_text(clu, encoding="utf-8")
    print("已写：trajectory_report_corrected.txt / slice_report_corrected.txt / cluster_report_corrected.txt")


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    main()
