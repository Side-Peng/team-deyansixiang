# -*- coding: utf-8 -*-
"""cluster_pm.py — 第四轮（扩展版）：PET−/Plasma+ 组内双亚群假说验证

合并症 = 双来源 ∪：MEDHIST 系统标记（ADNI1/GO/2）+ INITHEALTH 关键词匹配（ADNI3+）。
设计：聚类输入只用 血浆标志物 + 合并症（[pT217, AB42/40, CARDIO, ENDO]），
认知与轨迹只作验证变量（避免循环）。方法：KMeans k=2（标准化），k=3 敏感性。
验证：两簇的轨迹、基线认知、诊断构成、年龄是否可区分。

用法：python scripts/cluster_pm.py
输出：data/processed/cluster_report_v2.txt
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2_contingency
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

import config as C
from trajectory import per_patient_change

FEATS = ["PTAU217_bl", "AB_RATIO_bl", "CARDIO", "ENDO"]

INITHEALTH_KEYWORDS = {  # 列名 → (命中即阳性, 需排除的否定表述)
    "CARDIO": (["hypert", "high blood pressure", "hbp"], ["without diagnosis of hypertension"]),
    "ENDO": (["diabet"], []),
    "CKD": (["chronic kidney", "renal", "kidney disease", "ckd"], ["without renal involvement"]),
}


def make_comorbidity(sub: pd.DataFrame) -> pd.DataFrame:
    """合并症双来源：MEDHIST 系统标记 ∪ INITHEALTH 关键词（每 RID 任一命中即阳性）。"""
    mh = pd.read_csv(C.RAW_DIR / C.FILES["medhist"], low_memory=False)
    mh = mh.sort_values(["RID", "VISDATE"]).groupby("RID", as_index=False).first()
    com = pd.DataFrame({"RID": mh["RID"],
                        "CARDIO_MH": C.normalize_yes_no(mh["MH4CARD"]).fillna(0).astype(int),
                        "ENDO_MH": C.normalize_yes_no(mh["MH9ENDO"]).fillna(0).astype(int)})

    ih = pd.read_csv(C.RAW_DIR / C.FILES["inithealth"], low_memory=False)
    # fillna("") 后 astype(str)：兼容 pandas>=3（3.x 下 astype(str) 对 NaN 保留 NaN 而非 "nan"，
    # 导致后续 `k in s` 对 float 报错；fillna("") 与旧版 "nan" 均不命中关键词，结果等价）。
    ih["_low"] = (ih["IHSYMPTOM"].fillna("").astype(str).str.lower() + " " +
                  ih["IHDESC"].fillna("").astype(str).str.lower())
    for col, (kws, excls) in INITHEALTH_KEYWORDS.items():
        hit = ih["_low"].apply(lambda s: any(k in s for k in kws) and not any(e in s for e in excls))
        com[col + "_IH"] = ih.loc[hit].groupby("RID").size().reindex(com["RID"]).fillna(0).astype(int).values

    sub = sub.merge(com, on="RID", how="left")
    sub["CARDIO"] = ((sub["CARDIO_MH"].fillna(0) + sub["CARDIO_IH"].fillna(0)) > 0).astype(int)
    sub["ENDO"] = ((sub["ENDO_MH"].fillna(0) + sub["ENDO_IH"].fillna(0)) > 0).astype(int)
    sub["CKD"] = (sub["CKD_IH"].fillna(0) > 0).astype(int)  # CKD_IH 是命中行数，需二值化
    sub["COM_SRC"] = np.where(sub["CARDIO_MH"].fillna(0) + sub["ENDO_MH"].fillna(0) > 0, "MEDHIST",
                              np.where(sub["CARDIO_IH"].fillna(0) + sub["ENDO_IH"].fillna(0) > 0, "INITHEALTH", "无"))
    return sub


def main():
    sub = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    sub = make_comorbidity(sub)

    cdr = pd.read_csv(C.RAW_DIR / C.FILES["cdr"], low_memory=False)
    adas = pd.read_csv(C.RAW_DIR / C.FILES["adas"], low_memory=False)
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")
    sub = sub.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    sub["D_ADAS13_yr"] = sub["D_ADAS13"] / sub["YRS_ADAS"]

    pm = sub[(sub.PET_STATUS == 0) & (sub.PLASMA_STATUS == 1)].copy()
    pm = pm.dropna(subset=FEATS)
    lines = [f"PET−/Plasma+ 组：{len(pm)} 人（聚类输入齐全；纵向 {pm.D_ADAS13_yr.notna().sum()} 人）",
             f"合并症来源构成: MEDHIST={ (pm.COM_SRC=='MEDHIST').sum() } INITHEALTH={ (pm.COM_SRC=='INITHEALTH').sum() } 无={(pm.COM_SRC=='无').sum()}",
             f"合并症阳性率（双来源）: 心血管 {pm.CARDIO.mean():.1%} 内分泌 {pm.ENDO.mean():.1%} CKD {pm.CKD.mean():.1%}"]

    X = StandardScaler().fit_transform(pm[FEATS])

    # 主分析 k=2
    km2 = KMeans(n_clusters=2, n_init=10, random_state=42).fit(X)
    pm["CLUSTER2"] = km2.labels_
    sil2 = silhouette_score(X, km2.labels_)
    lines.append(f"\nKMeans k=2: 轮廓系数={sil2:.3f}  簇规模: {np.bincount(km2.labels_).tolist()}")

    # 两簇特征与验证变量对比
    for cl in sorted(pm.CLUSTER2.unique()):
        d = pm[pm.CLUSTER2 == cl]
        tr = d.D_ADAS13_yr.dropna()
        lines.append(f"\n[簇 {cl}] n={len(d)}")
        for f in FEATS + ["AGE", "ADAS13_bl"]:
            if f in d:
                lines.append(f"  {f:<12} 中位={d[f].median():.3f}")
        lines.append(f"  心血管阳性 {d.CARDIO.mean():.1%} / 内分泌 {d.ENDO.mean():.1%} / CKD {d.CKD.mean():.1%}")
        lines.append(f"  ΔADAS13/年 中位={tr.median():+.2f} (n={len(tr)})")
        lines.append(f"  诊断构成: CN { (d.DX_bl==1).mean():.0%} MCI {(d.DX_bl==2).mean():.0%} AD {(d.DX_bl==3).mean():.0%}")
        lines.append(f"  年龄中位 {d.AGE.median():.1f}")

    # 两簇统计检验
    c0 = pm[pm.CLUSTER2 == 0]; c1 = pm[pm.CLUSTER2 == 1]
    lines.append("\n===== 两簇判别检验 =====")
    lines.append(f"轨迹 ΔADAS13/年: 簇0={c0.D_ADAS13_yr.median():+.2f} vs 簇1={c1.D_ADAS13_yr.median():+.2f}, "
                 f"Mann-Whitney p={mannwhitneyu(c0.D_ADAS13_yr.dropna(), c1.D_ADAS13_yr.dropna()).pvalue:.4g}")
    lines.append(f"基线 ADAS13: 簇0={c0.ADAS13_bl.median():.1f} vs 簇1={c1.ADAS13_bl.median():.1f}, "
                 f"p={mannwhitneyu(c0.ADAS13_bl.dropna(), c1.ADAS13_bl.dropna()).pvalue:.4g}")
    for col, name in [("CARDIO", "心血管"), ("ENDO", "内分泌")]:
        tab = [[(c0[col] == 1).sum(), (c0[col] == 0).sum()],
               [(c1[col] == 1).sum(), (c1[col] == 0).sum()]]
        chi2, p, _, _ = chi2_contingency(np.array(tab))
        lines.append(f"{name}构成: 簇0={c0[col].mean():.1%} 簇1={c1[col].mean():.1%}, 卡方 p={p:.4g}")
    tab2 = [[(c0["CKD"] == 1).sum(), (c0["CKD"] == 0).sum()],
            [(c1["CKD"] == 1).sum(), (c1["CKD"] == 0).sum()]]
    chi2, p_ckd, _, _ = chi2_contingency(np.array(tab2))
    lines.append(f"CKD构成: 簇0={c0.CKD.mean():.1%} 簇1={c1.CKD.mean():.1%}, 卡方 p={p_ckd:.4g}")

    # k=3 敏感性
    km3 = KMeans(n_clusters=3, n_init=10, random_state=42).fit(X)
    sil3 = silhouette_score(X, km3.labels_)
    pm["CLUSTER3"] = km3.labels_
    lines.append(f"\nKMeans k=3 敏感性: 轮廓系数={sil3:.3f} 簇规模={np.bincount(km3.labels_).tolist()}")
    for cl in sorted(pm.CLUSTER3.unique()):
        d = pm[pm.CLUSTER3 == cl]
        lines.append(f"  [簇 {cl}] n={len(d)} 合并症率 心/内/CKD={d.CARDIO.mean():.0%}/{d.ENDO.mean():.0%}/{d.CKD.mean():.0%} "
                     f"ΔADAS13/年={d.D_ADAS13_yr.median():+.2f} 基线ADAS13={d.ADAS13_bl.median():.1f}")

    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "cluster_report_v2.txt").write_text(txt, encoding="utf-8")
    print("\n已保存: data/processed/cluster_report_v2.txt（对比 v1：合并症双来源，样本应显著扩大）")


if __name__ == "__main__":
    main()
