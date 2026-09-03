# -*- coding: utf-8 -*-
"""slice_analysis.py — 第三轮：合并症混杂检查 + 诊断分层复现性

① 合并症代理（MEDHIST 系统标记：MH4CARD 心血管、MH9ENDO 内分泌）
   - 四组阳性率 + 卡方检验
   - 轨迹回归加入合并症后，PET−/Plasma+ 系数是否保留
② 诊断分层（CN/MCI/AD）：层内重跑四组结构，检验"中间态"是否跨层复现

用法：python scripts/slice_analysis.py
输出：data/processed/slice_report.txt
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu

import config as C

DX_LABEL = {1.0: "CN", 2.0: "MCI", 3.0: "AD", 10.0: "UNKNOWN(10)"}


def main():
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    sub = sub.dropna(subset=["PET_STATUS", "PLASMA_STATUS", "AGE"]).copy()
    sub["GROUP"] = (sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int))
    labels = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}

    # 合并症代理：每 RID 最早访问
    mh = pd.read_csv(C.RAW_DIR / C.FILES["medhist"], low_memory=False)
    mh[C.FIELDS["date"]["ptdemog"]] = pd.to_datetime(mh["VISDATE"], errors="coerce")
    mh = mh.sort_values(["RID", "VISDATE"]).groupby("RID", as_index=False).first()
    sub = sub.merge(mh[["RID", "MH4CARD", "MH9ENDO"]], on="RID", how="left")
    sub["CARDIO"] = C.normalize_yes_no(sub["MH4CARD"]).fillna(0).astype(int)
    sub["ENDO"] = C.normalize_yes_no(sub["MH9ENDO"]).fillna(0).astype(int)

    lines = [f"样本：{len(sub)}（MEDHIST 覆盖 {sub.MH4CARD.notna().sum()} 人，ADNI1/GO/2）"]

    # ---------- ① 合并症 ----------
    lines.append("\n===== ① 合并症代理（MEDHIST 系统标记）=====")
    for com, name in [("CARDIO", "心血管 MH4CARD"), ("ENDO", "内分泌 MH9ENDO")]:
        lines.append(f"\n{name} 阳性率：")
        tbl = []
        for g in sorted(sub.GROUP.unique()):
            d = sub[sub.GROUP == g]
            r = d[com].mean()
            tbl.append([d[com].sum(), len(d) - d[com].sum()])
            lines.append(f"  {labels[g]:<14} {r:.2%}  (n={len(d)})")
        chi2, p, _, _ = chi2_contingency(np.array(tbl))
        lines.append(f"  卡方检验: chi2={chi2:.1f}, p={p:.4g}")

    # 轨迹回归加合并症（复用 trajectory 数据构建）
    cdr = pd.read_csv(C.RAW_DIR / C.FILES["cdr"], low_memory=False)
    adas = pd.read_csv(C.RAW_DIR / C.FILES["adas"], low_memory=False)
    from trajectory import per_patient_change
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")
    m = sub.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    m["D_ADAS13_yr"] = m["D_ADAS13"] / m["YRS_ADAS"]

    from sklearn.linear_model import LinearRegression
    for col in ["D_ADAS13_yr"]:
        d = m.dropna(subset=[col, "AGE", "GENDER", "EDUCAT", "MH4CARD"])
        d["G1"] = (d.GROUP == 1).astype(int); d["G2"] = (d.GROUP == 2).astype(int); d["G3"] = (d.GROUP == 3).astype(int)
        X = d[["G1", "G2", "G3", "AGE", "GENDER", "EDUCAT", "CARDIO", "ENDO"]].copy()
        X["GENDER"] = X["GENDER"].map({1: 1, 2: 0})
        X["CARDIO"] = X["CARDIO"].astype(int); X["ENDO"] = X["ENDO"].astype(int)
        lm = LinearRegression().fit(X, d.loc[X.index, col])
        lines.append(f"\n{col}: 校正年龄/性别/教育 + 合并症后组系数（相对 PET−/Plasma−）")
        for nm, coef in zip(["PET−/Plasma+", "PET+/Plasma−", "PET+/Plasma+", "CARDIO", "ENDO"],
                            list(lm.coef_[:3]) + [lm.coef_[6], lm.coef_[7]]):
            lines.append(f"  {nm:<14}: {coef:+.3f}")

    # ---------- ② 诊断分层 ----------
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
            adas_bl = d.ADAS13_bl.median()
            traj = d.D_ADAS13_yr.median()
            row.append((labels[g], len(d), adas_bl, traj))
        lines.append(f"  {'组':<14} {'n':>4} {'ADAS13_bl':>9} {'ΔADAS13/年':>9}")
        for g, n, ab, tr in row:
            lines.append(f"  {g:<14} {n:>4} {ab:>9.1f} {tr:>+9.2f}")
        # 中间态检验：层内 PET−/Plasma+（GROUP==1）vs 双阴（轨迹）
        p2 = layer[layer.GROUP == 1].D_ADAS13_yr.dropna().values
        p0 = layer[layer.GROUP == 0].D_ADAS13_yr.dropna().values
        p3 = layer[layer.GROUP == 3].D_ADAS13_yr.dropna().values
        if len(p2) >= 10 and len(p0) >= 10 and len(p3) >= 10:
            p_vs_neg = mannwhitneyu(p2, p0).pvalue
            p_vs_pos = mannwhitneyu(p2, p3).pvalue
            med_neg, med_2, med_pos = np.median(p0), np.median(p2), np.median(p3)
            strict_mid = med_neg <= med_2 <= med_pos and p_vs_neg < 0.05
            lines.append(f"  中间态判定: 中位 双阴={med_neg:+.2f} ≤ P−/P+={med_2:+.2f} ≤ 双阳={med_pos:+.2f}")
            lines.append(f"  P−/P+ vs 双阴 p={p_vs_neg:.4g}；vs 双阳 p={p_vs_pos:.4g} → {'✅ 严格中间态' if strict_mid else '❌ 非严格中间态'}")

    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "slice_report.txt").write_text(txt, encoding="utf-8")
    print("\n已保存: data/processed/slice_report.txt")


if __name__ == "__main__":
    main()
