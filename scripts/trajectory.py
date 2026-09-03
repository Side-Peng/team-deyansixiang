# -*- coding: utf-8 -*-
"""trajectory.py — 第二轮：四组纵向认知轨迹比较

基线=各表最早访问，随访终点=最晚访问；Δ 年化 = (last−first)/随访年数。
步骤：① 年化变化率统计；② Kruskal-Wallis + 两两 Mann-Whitney；
③ 年龄/性别/教育校正回归（sklearn LinearRegression）；
④ 置换检验（打乱组标签 1000 次，KW H 统计量）→ 排除随机运气。

用法：python scripts/trajectory.py
输出：data/processed/trajectory_report.txt
"""
import numpy as np
import pandas as pd
import config as C


def per_patient_change(df: pd.DataFrame, date_col: str, val_col: str, out_delta: str, out_years: str) -> pd.DataFrame:
    """每 RID：最早与最晚访问的 Δ 及随访年数。"""
    d = df.copy()
    d[date_col] = pd.to_datetime(d[date_col], errors="coerce")
    d[val_col] = pd.to_numeric(d[val_col], errors="coerce")
    d = d.dropna(subset=[date_col, val_col])
    d = d.sort_values(["RID", date_col])
    g = d.groupby("RID", as_index=False)[[date_col, val_col]]
    first = g.first(); last = g.last()
    years = (last[date_col] - first[date_col]).dt.days / 365.25
    return pd.DataFrame({
        "RID": first["RID"],
        out_delta: last[val_col].values - first[val_col].values,
        out_years: years.values,
    })


def kruskal_h(y, group):
    from scipy.stats import kruskal, mannwhitneyu
    groups = [y[group == g] for g in np.unique(group)]
    H, p = kruskal(*groups)
    return H, p, mannwhitneyu


def main():
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    sub = sub.dropna(subset=["PET_STATUS", "PLASMA_STATUS", "AGE"]).copy()
    sub["GROUP"] = (sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int))
    labels = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}

    cdr = pd.read_csv(C.RAW_DIR / C.FILES["cdr"], low_memory=False)
    adas = pd.read_csv(C.RAW_DIR / C.FILES["adas"], low_memory=False)
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")

    m = sub.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    m["D_CDRSB_yr"] = m["D_CDRSB"] / m["YRS_CDR"]
    m["D_ADAS13_yr"] = m["D_ADAS13"] / m["YRS_ADAS"]

    from scipy.stats import kruskal, mannwhitneyu
    lines = [f"样本：{len(m)}（基线双状态+年龄齐全）；随访: CDRSB {m.YRS_CDR.notna().sum()} 人 / ADAS13 {m.YRS_ADAS.notna().sum()} 人"]

    # ① 年化变化率
    lines.append("\n===== ① 年化变化率（Δ/年，中位数 [IQR]）=====")
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        lines.append(f"\n{col}:")
        for g in sorted(m.GROUP.unique()):
            v = m.loc[m.GROUP == g, col].dropna()
            q = v.quantile([.5, .25, .75])
            lines.append(f"  {labels[g]:<14} n={len(v):>4}  中位={q[.5]:+.2f} [{q[.25]:+.2f}, {q[.75]:+.2f}]  随访年数={m.loc[m.GROUP==g,'YRS_CDR'].median():.1f}")

    # ② KW + 两两（重点 PET−/Plasma+ vs 双阴、双阳）
    lines.append("\n===== ② 组间检验（Kruskal-Wallis + Mann-Whitney）=====")
    pairs = [(1, 0, "PET−/Plasma+ vs PET−/Plasma−"),
             (1, 3, "PET−/Plasma+ vs PET+/Plasma+"),
             (1, 2, "PET−/Plasma+ vs PET+/Plasma−")]
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        d = m.dropna(subset=[col])
        H, p = kruskal(*[d.loc[d.GROUP == g, col].values for g in sorted(d.GROUP.unique())])
        lines.append(f"\n{col}: KW H={H:.1f}, p={p:.4g}")
        for g1, g2, name in pairs:
            a1 = d.loc[d.GROUP == g1, col].values; a2 = d.loc[d.GROUP == g2, col].values
            if len(a1) and len(a2):
                p2 = mannwhitneyu(a1, a2, alternative="two-sided").pvalue
                lines.append(f"  {name}: n={len(a1)}/{len(a2)} p={p2:.4g}")

    # ③ 混杂校正（Δ/年 ~ group + age + sex + educ）
    lines.append("\n===== ③ 混杂校正（LinearRegression, 以 PET−/Plasma− 为参照）=====")
    from sklearn.linear_model import LinearRegression
    m["G1"] = (m.GROUP == 1).astype(int)  # 真 PET−/Plasma+
    m["G2"] = (m.GROUP == 2).astype(int)  # 真 PET+/Plasma−
    m["G3"] = (m.GROUP == 3).astype(int)
    for col in ["D_CDRSB_yr", "D_ADAS13_yr"]:
        d = m.dropna(subset=[col, "AGE", "GENDER", "EDUCAT"])
        X = d[["G1", "G2", "G3", "AGE", "GENDER", "EDUCAT"]].copy()
        X["GENDER"] = X["GENDER"].map({1: 1, 2: 0})  # ADNI 数字编码: 1=Male, 2=Female
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
        y = d[col].values; g = d.GROUP.values
        H_obs = kruskal(*[y[g == x] for x in np.unique(g)])[0]
        cnt = 0
        for _ in range(1000):
            gp = rng.permutation(g)
            H_p = kruskal(*[y[gp == x] for x in np.unique(g)])[0]
            cnt += H_p >= H_obs
        lines.append(f"{col}: 观察 H={H_obs:.1f}，置换 H≥观察 比例 p={cnt/1000:.3f}")

    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "trajectory_report.txt").write_text(txt, encoding="utf-8")
    print("\n已保存: data/processed/trajectory_report.txt")


if __name__ == "__main__":
    main()
