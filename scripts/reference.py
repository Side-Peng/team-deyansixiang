# -*- coding: utf-8 -*-
"""reference.py — 参照系设计（GOAI AI4S 探索赛道 · 复赛三件套之二）

对齐 four_page_guide.md Page 3 §3.2 的四层参照，给每个「发现」一个可运行的
对照框架，回答「这个发现比平凡/随机/不干预/Pyun baseline 多了什么」。

四层参照：
  ① 平凡解 (trivial)         —— 无切分（全样本单一结构）；单模态二分（仅血浆 / 仅影像）
  ② 随机参照 (random)        —— 置换检验：真实分离度须显著超出置换分布
  ③ 无干预参照 (no-intervention)—— 不定义不切分的原始四组结构（未做亚型发现前）
  ④ 非平凡 baseline (nontrivial)—— Pyun 2024 四组结构；发现须落在「四组内部亚结构 /
                                     跨切片复现性 / 混杂归因」三层中至少一层

用法：
    from environment import ExplorationEnvironment
    from reference import ReferenceFrame
    env = ExplorationEnvironment()
    ref = ReferenceFrame(env)
    report = ref.report()                      # → dict
    ref.write_report()                          # → data/processed/reference_report.txt
"""
import numpy as np
import pandas as pd
from scipy import stats as sps

import config as C
import environment as E
from environment import (compute_plasma_status, add_group, build_longitudinal,
                         GROUP_LABELS, OUTCOME_COL, _cohens_d_ci, _permutation_p,
                         _cross_slice_reproducibility, _adjusted_coeffs, _clean)


class ReferenceFrame:
    def __init__(self, env: "E.ExplorationEnvironment", assay_pair="ptau217_vs_pet",
                 missing_policy="exclude", outcome="ADAS13"):
        self.env = env
        self.assay_pair = assay_pair
        self.missing_policy = missing_policy
        self.outcome = OUTCOME_COL[outcome]
        self._frame = self._build_frame()

    def _build_frame(self) -> pd.DataFrame:
        thr = (C.THRESHOLDS["p_tau217_pg_ml"] if self.assay_pair == "ptau217_vs_pet"
               else C.THRESHOLDS["abeta_ratio"])
        plasma = compute_plasma_status(self.env.wide, self.assay_pair, thr, self.missing_policy)
        sub = add_group(self.env.wide, plasma)
        m = build_longitudinal(sub)
        from cluster_pm import make_comorbidity
        m = make_comorbidity(m)
        m["GENDER"] = m["GENDER"].map({1: 1, 2: 0})
        return m

    # ---- ① 平凡解 ----
    def trivial(self) -> dict:
        m = self._frame.dropna(subset=[self.outcome])
        whole = m[self.outcome]
        # 单模态二分：仅 PET（PET+ vs PET−）、仅血浆（Plasma+ vs Plasma−）
        pet_pos = m[m.PET_STATUS == 1][self.outcome]
        pet_neg = m[m.PET_STATUS == 0][self.outcome]
        plasma_pos = m[m.PLASMA_STATUS == 1][self.outcome]
        plasma_neg = m[m.PLASMA_STATUS == 0][self.outcome]
        def _med_iqr(v):
            v = v.dropna()
            q = v.quantile([.5, .25, .75])
            return {"median": float(q[.5]), "iqr": [float(q[.25]), float(q[.75])], "n": int(len(v))}
        return {
            "whole_sample_no_split": _med_iqr(whole),
            "pet_only_binary": {"pet_pos": _med_iqr(pet_pos), "pet_neg": _med_iqr(pet_neg),
                                 "cohens_d": _cohens_d_ci(pet_pos, pet_neg),
                                 "mw_p": float(sps.mannwhitneyu(pet_pos.dropna(), pet_neg.dropna()).pvalue)},
            "plasma_only_binary": {"plasma_pos": _med_iqr(plasma_pos), "plasma_neg": _med_iqr(plasma_neg),
                                    "cohens_d": _cohens_d_ci(plasma_pos, plasma_neg),
                                    "mw_p": float(sps.mannwhitneyu(plasma_pos.dropna(), plasma_neg.dropna()).pvalue)},
        }

    # ---- ② 随机参照 ----
    def random(self, n_perm=500) -> dict:
        m = self._frame.dropna(subset=["GROUP", self.outcome])
        H_obs, p_perm = _permutation_p(m, "GROUP", self.outcome, n_perm=n_perm, seed=self.env.seed)
        # 置换分布分位数（重算一次拿分位数，避免与 p 值重复跑时口径不一致——这里直接给出）
        d = m.dropna(subset=["GROUP", self.outcome])
        y = d[self.outcome].values
        g = d.GROUP.values
        uniq = np.unique(g)
        rng = np.random.default_rng(self.env.seed)
        Hperm = []
        for _ in range(n_perm):
            gp = rng.permutation(g)
            Hperm.append(sps.kruskal(*[y[gp == x] for x in uniq])[0])
        return {"H_obs": H_obs, "perm_p": p_perm,
                "perm_dist": {"p95": float(np.percentile(Hperm, 95)),
                              "p99": float(np.percentile(Hperm, 99)),
                              "max": float(np.max(Hperm))},
                "verdict": "超出随机参照" if (H_obs is not None and H_obs > np.percentile(Hperm, 95)) else "未超随机参照"}

    # ---- ③ 无干预参照（原始四组，未做亚型发现） ----
    def no_intervention(self) -> dict:
        m = self._frame.dropna(subset=["GROUP", self.outcome])
        med = {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, self.outcome].median()) for g in GROUP_LABELS}
        p_int_neg, _, _ = E._mannwhitney(m, "GROUP", self.outcome, 1, 0)
        p_int_pos, _, _ = E._mannwhitney(m, "GROUP", self.outcome, 1, 3)
        return {"four_group_medians": med,
                "intermediate_vs_neg": {"p": p_int_neg,
                                         "cohens_d": _cohens_d_ci(m.loc[m.GROUP == 1, self.outcome],
                                                                   m.loc[m.GROUP == 0, self.outcome])},
                "intermediate_vs_pos": {"p": p_int_pos}}

    # ---- ④ 非平凡 baseline（Pyun 四组）+ 三层增值判定 ----
    def nontrivial(self) -> dict:
        m = self._frame
        # 跨切片复现性（增值层 2）
        repro = _cross_slice_reproducibility(m, self.outcome)
        # 混杂归因（增值层 3）：加入合并症后中间态系数是否保留
        base = _adjusted_coeffs(m, "GROUP", self.outcome, ["AGE", "GENDER", "EDUCAT"], ref_val=0)
        full = _adjusted_coeffs(m, "GROUP", self.outcome, ["AGE", "GENDER", "EDUCAT", "CARDIO", "ENDO", "CKD"], ref_val=0)
        base_g1 = base["group_coefs"][1] if base["group_coefs"] else None
        full_g1 = full["group_coefs"][1] if full["group_coefs"] else None
        delta = (full_g1 - base_g1) if (base_g1 is not None and full_g1 is not None) else None
        # 增值层 1（四组内部亚结构）在 discover_subtypes 动作里体现，这里给出判定口径说明
        return {
            "pyun_four_group_baseline": self.no_intervention(),
            "add_value_layer_2_cross_slice": {"repro_coeff": repro["coeff"],
                                               "per_slice": repro["per_slice"]},
            "add_value_layer_3_confounder": {"base_coef_G1": base_g1,
                                              "comorbidity_adjusted_coef_G1": full_g1,
                                              "delta": delta},
            "add_value_layer_1_within_group_substructure": {"note": "由 discover_subtypes 动作评估（组内聚类是否给出轨迹可区分的亚型）"},
        }

    def report(self) -> dict:
        return _clean({
            "assay_pair": self.assay_pair,
            "missing_policy": self.missing_policy,
            "outcome": self.outcome,
            "trivial": self.trivial(),
            "random": self.random(),
            "no_intervention": self.no_intervention(),
            "nontrivial": self.nontrivial(),
        })

    def write_report(self) -> "Path":
        import json
        rep = self.report()
        out = C.PROC_DIR / "reference_report.json"
        out.write_text(json.dumps(rep, ensure_ascii=False, indent=1), encoding="utf-8")
        txt = self._format_text(rep)
        (C.PROC_DIR / "reference_report.txt").write_text(txt, encoding="utf-8")
        return C.PROC_DIR / "reference_report.txt"

    def _format_text(self, rep: dict) -> str:
        lines = ["参照系报告（Page 3 §3.2 四层参照）",
                 f"口径：{rep['assay_pair']} / missing_policy={rep['missing_policy']} / outcome={rep['outcome']}", ""]
        t = rep["trivial"]
        lines.append("① 平凡解：")
        lines.append(f"   无切分全样本：{t['whole_sample_no_split']}")
        lines.append(f"   单模态-仅PET：d={t['pet_only_binary']['cohens_d']['d']:.3f} p={t['pet_only_binary']['mw_p']:.4g}")
        lines.append(f"   单模态-仅血浆：d={t['plasma_only_binary']['cohens_d']['d']:.3f} p={t['plasma_only_binary']['mw_p']:.4g}")
        r = rep["random"]
        lines.append(f"② 随机参照：H_obs={r['H_obs']:.1f}，置换 p95={r['perm_dist']['p95']:.1f}，perm_p={r['perm_p']:.4g} → {r['verdict']}")
        n = rep["no_intervention"]
        lines.append(f"③ 无干预参照（原始四组）：{n['four_group_medians']}")
        nt = rep["nontrivial"]
        lines.append(f"④ 非平凡 baseline（Pyun 四组）：中间态 vs 双阴 p={nt['pyun_four_group_baseline']['intermediate_vs_neg']['p']:.4g}")
        lines.append(f"   增值层2-跨切片复现系数={nt['add_value_layer_2_cross_slice']['repro_coeff']}")
        lines.append(f"   增值层3-混杂归因：base_G1={nt['add_value_layer_3_confounder']['base_coef_G1']} "
                     f"校正后_G1={nt['add_value_layer_3_confounder']['comorbidity_adjusted_coef_G1']} "
                     f"Δ={nt['add_value_layer_3_confounder']['delta']}")
        return "\n".join(lines)


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    env = E.ExplorationEnvironment()
    ref = ReferenceFrame(env)
    print(ref.write_report().read_text(encoding="utf-8"))
