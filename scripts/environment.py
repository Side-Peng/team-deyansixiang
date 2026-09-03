# -*- coding: utf-8 -*-
"""environment.py — 最小可运行 Agent 探索环境（GOAI AI4S 探索赛道 · 复赛核心件）

把初赛阶段的静态分析管线（preprocess → calibrate → define_discordance →
trajectory → slice_analysis → cluster_pm → tau_analysis）封装为一个统一的
「观察 / 行动 / 反馈 / 记录」接口，供 Agent 进入、尝试、受挫、修正并持续探索。

接口（与 four_page_guide.md Page 2 §2.2 对应）：
    env = ExplorationEnvironment()
    obs  = env.observe()                 # 观察：可读字段 + 当前状态
    fb   = env.act("define_discordance", {"assay_pair": "ptau217_vs_pet"})
                                         # 行动：5 个动作之一，返回 5 项反馈
    env.record(round_id, action, params, rationale, fb, state_update)
                                         # 记录：round_log 落盘 data/processed/exploration_log.jsonl

五个动作（Page 2 §2.2 行动）：
    define_discordance(assay_pair, plasma_threshold, outcome)
    discover_subtypes(method, k, features)
    select_slice(diagnosis, age_band)
    test_confounder(control_vars)
    sensitivity_analysis(outcome, window)

每轮反馈统一返回 5 项（Page 2 §2.2 反馈）：
    ① discordant_group_size      不一致组规模与构成
    ② cross_slice_reproducibility  跨切片复现性系数
    ③ trajectory_separation       亚型/组间认知轨迹差异（效应量 + CI + p）
    ④ confounder_adjusted_change  混杂控制后效应变化
    ⑤ vs_random_reference         与随机参照的差异量

说明：本环境只依赖 data/processed/subjects_wide.csv（初赛管线已生成）+ 原始表；
不重写初赛脚本，只在其上做一层统一接口，保证可逐轮复现。
"""
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score

import config as C
from trajectory import per_patient_change
from cluster_pm import make_comorbidity
from data_layer import ALIGNMENT_WINDOWS

# 四组命名（GROUP = PET*2 + PLASMA，与初赛约定一致，见 CLAUDE.md 关键约定）
GROUP_LABELS = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}

DX_LABEL = {1.0: "CN", 2.0: "MCI", 3.0: "AD", 10.0: "UNKNOWN(10)"}

OUTCOME_COL = {"ADAS13": "D_ADAS13_yr", "CDRSB": "D_CDRSB_yr"}

# 动作参数空间（供 Agent 观察 + 校验参数合法性）
ACTION_SCHEMAS = {
    "define_discordance": {
        "assay_pair": ["ptau217_vs_pet", "abeta_ratio_vs_pet"],
        "plasma_threshold": "float（可选，覆盖 config 默认阈值）",
        "outcome": ["ADAS13", "CDRSB"],
        "missing_policy": ["exclude", "missing_as_negative"],
    },
    "discover_subtypes": {
        "method": ["kmeans", "gmm"],
        "k": [2, 3, 4],
        "diagnosis": ["ALL", "CN", "MCI", "AD"],
        "features": ["PTAU217_bl", "AB_RATIO_bl", "CARDIO", "ENDO", "CKD"],
    },
    "select_slice": {
        "diagnosis": ["CN", "MCI", "AD", "ALL"],
        "age_band": ["<75", ">=75", "ALL"],
    },
    "test_confounder": {
        "control_vars": ["AGE", "GENDER", "EDUCAT", "CARDIO", "ENDO", "CKD"],
    },
    "sensitivity_analysis": {
        "outcome": ["ADAS13", "CDRSB"],
        "window": ["all", ">=2yr_followup"],
    },
    "profile_mechanism": {
        "target_group": ["PET−/Plasma+", "PET+/Plasma−"],
        "markers": ["tau_pet", "GFAP", "NfL"],
        "alignment_window": ["90d", "180d", "365d"],
        "outcome": ["ADAS13", "CDRSB"],
    },
}

# 探索日志落盘路径
LOG_PATH = C.PROC_DIR / "exploration_log.jsonl"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _clean(v):
    """递归把 numpy 标量转成 Python 原生类型，保证 JSON 可序列化。"""
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.ndarray,)):
        return _clean(v.tolist())
    return v


# --------------------------------------------------------------------------
# 数据装载与派生
# --------------------------------------------------------------------------
def load_wide() -> pd.DataFrame:
    p = C.PROC_DIR / "subjects_wide.csv"
    if not p.exists():
        raise SystemExit("[environment] 缺少 subjects_wide.csv：先跑 preprocess.py + calibrate.py")
    return pd.read_csv(p)


def compute_plasma_status(df: pd.DataFrame, assay_pair: str, threshold: float,
                          missing_policy: str = "exclude") -> pd.Series:
    """按标志物对 + 阈值计算血浆阳性状态。逆向标志物（Aβ42/40 比值越低越病理）取 < 阈值。

    missing_policy:
      "exclude"             —— 缺失标志物 → NaN（正确口径：没测 ≠ 阴性；由 add_group 排除）
      "missing_as_negative" —— 缺失标志物 → 0（初赛 preprocess.py 的历史口径，会系统性
                               把「未测血浆」混进「血浆阴性」，用于探索环节复现/证伪该 bug）
    """
    if assay_pair == "ptau217_vs_pet":
        col = "PTAU217_bl"
        pos = df[col] > threshold
    elif assay_pair == "abeta_ratio_vs_pet":
        col = "AB_RATIO_bl"
        pos = df[col] < threshold
    else:
        raise ValueError(f"未知 assay_pair: {assay_pair}")
    status = pos.astype("Int64")
    status[df[col].isna()] = pd.NA
    if missing_policy == "missing_as_negative":
        status = status.fillna(0)
    return status


def add_group(df: pd.DataFrame, plasma_status: pd.Series) -> pd.DataFrame:
    d = df.copy()
    d = d.dropna(subset=["PET_STATUS"]).copy()
    d["PLASMA_STATUS"] = plasma_status.loc[d.index]
    d = d.dropna(subset=["PET_STATUS", "PLASMA_STATUS"]).copy()
    d["GROUP"] = (d["PET_STATUS"].astype(int) * 2 + d["PLASMA_STATUS"].astype(int))
    return d


def build_longitudinal(sub: pd.DataFrame) -> pd.DataFrame:
    """宽表 + 纵向年化变化率（ADAS13 / CDRSB）。复用初赛 trajectory.per_patient_change。"""
    cdr = pd.read_csv(C.RAW_DIR / C.FILES["cdr"], low_memory=False)
    adas = pd.read_csv(C.RAW_DIR / C.FILES["adas"], low_memory=False)
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")
    m = sub.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    m["D_CDRSB_yr"] = m["D_CDRSB"] / m["YRS_CDR"]
    m["D_ADAS13_yr"] = m["D_ADAS13"] / m["YRS_ADAS"]
    return m


# --------------------------------------------------------------------------
# 统计原语（每轮反馈用）
# --------------------------------------------------------------------------
def _cohens_d_ci(a, b, n_boot=300, seed=2026):
    """Cohen's d（合并标准差）+ 自助法 95% CI。"""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return {"d": None, "ci_lo": None, "ci_hi": None, "n_a": int(len(a)), "n_b": int(len(b))}
    na, nb = len(a), len(b)
    var_a, var_b = a.var(ddof=1), b.var(ddof=1)
    sp = np.sqrt(((na - 1) * var_a + (nb - 1) * var_b) / (na + nb - 2))
    d = (a.mean() - b.mean()) / sp if sp > 0 else 0.0
    rng = np.random.default_rng(seed)
    ds = []
    for _ in range(n_boot):
        aa = rng.choice(a, na, replace=True)
        bb = rng.choice(b, nb, replace=True)
        spp = np.sqrt(((na - 1) * aa.var(ddof=1) + (nb - 1) * bb.var(ddof=1)) / (na + nb - 2))
        ds.append((aa.mean() - bb.mean()) / spp if spp > 0 else 0.0)
    lo, hi = np.percentile(ds, [2.5, 97.5])
    return {"d": float(d), "ci_lo": float(lo), "ci_hi": float(hi), "n_a": int(na), "n_b": int(nb)}


def _kw(m, group_col, outcome):
    """Kruskal-Wallis 跨组检验。返回 (H, p, cleaned_df)。"""
    d = m.dropna(subset=[group_col, outcome]).copy()
    if d.empty:
        return None, None, d
    gs = [d.loc[d[group_col] == g, outcome].values for g in sorted(d[group_col].unique())]
    H, p = sps.kruskal(*gs)
    return float(H), float(p), d


def _mannwhitney(m, group_col, outcome, g1, g2):
    d = m.dropna(subset=[group_col, outcome])
    a = d.loc[d[group_col] == g1, outcome].values
    b = d.loc[d[group_col] == g2, outcome].values
    if len(a) == 0 or len(b) == 0:
        return None, len(a), len(b)
    p = sps.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return float(p), int(len(a)), int(len(b))


def _permutation_p(m, group_col, outcome, n_perm=500, seed=2026):
    """置换检验：打乱组标签重算 KW H，返回 (H_obs, p_perm)。"""
    d = m.dropna(subset=[group_col, outcome])
    if d.empty:
        return None, None
    y = d[outcome].values
    g = d[group_col].values
    uniq = np.unique(g)
    if len(uniq) < 2:
        return None, None
    H_obs = sps.kruskal(*[y[g == x] for x in uniq])[0]
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        gp = rng.permutation(g)
        Hp = sps.kruskal(*[y[gp == x] for x in uniq])[0]
        cnt += Hp >= H_obs
    return float(H_obs), float((cnt + 1) / (n_perm + 1))


def _cluster_stability(X, labels, method, k, n_boot=25, seed=2026):
    if len(X) < max(20, k * 5):
        return {"n_boot": 0, "ari_mean": None, "ari_ci95": [None, None]}
    rng = np.random.default_rng(seed)
    scores = []
    for _ in range(n_boot):
        idx = np.sort(rng.choice(len(X), size=max(k * 5, int(len(X) * 0.8)), replace=False))
        if method == "kmeans":
            model = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X[idx])
            sampled = model.labels_
        else:
            model = GaussianMixture(n_components=k, random_state=seed).fit(X[idx])
            sampled = model.predict(X[idx])
        scores.append(float(adjusted_rand_score(labels[idx], sampled)))
    return {"n_boot": int(len(scores)), "ari_mean": float(np.mean(scores)),
            "ari_ci95": [float(x) for x in np.percentile(scores, [2.5, 97.5])]}


def _adjusted_coeffs(m, group_col, outcome, control_vars, ref_val=0):
    """线性回归：outcome ~ group(one-hot) + control_vars。返回组系数 dict + 完整模型 R²。"""
    d = m.dropna(subset=[group_col, outcome] + control_vars).copy()
    if d.empty:
        return {"group_coefs": None, "control_coefs": None, "r2": None}
    groups = sorted(int(g) for g in d[group_col].unique())
    # 组 one-hot 用字符串列名 "G{int}"，避免 sklearn>=1.9 对「int+str 混合列名」报错
    X = pd.DataFrame({f"G{int(g)}": (d[group_col].astype(int) == g).astype(int)
                      for g in groups if g != ref_val})
    for cv in control_vars:
        X[cv] = pd.to_numeric(d[cv], errors="coerce")
    X = X.dropna()
    if X.empty:
        return {"group_coefs": None, "control_coefs": None, "r2": None}
    y = d.loc[X.index, outcome]
    lm = LinearRegression().fit(X, y)
    names = list(X.columns)
    coefs = {int(ref_val): 0.0}
    for g in groups:
        if g != ref_val:
            coefs[int(g)] = float(lm.coef_[names.index(f"G{int(g)}")])
    ctrl_coefs = {cv: float(lm.coef_[names.index(cv)]) for cv in control_vars if cv in names}
    return {"group_coefs": coefs, "control_coefs": ctrl_coefs, "r2": float(lm.score(X, y))}


def _intermediate_check(m, outcome):
    """检验「双阴 ≤ 中间态 ≤ 双阳」的严格中间态是否成立（对齐 slice_analysis.py）。"""
    d = m.dropna(subset=["GROUP", outcome])
    p0 = d.loc[d.GROUP == 0, outcome].dropna()
    p1 = d.loc[d.GROUP == 1, outcome].dropna()
    p3 = d.loc[d.GROUP == 3, outcome].dropna()
    if len(p0) < 10 or len(p1) < 10 or len(p3) < 10:
        return {"holds": None, "n0": int(len(p0)), "n1": int(len(p1)), "n3": int(len(p3))}
    m0, m1, m3 = np.median(p0), np.median(p1), np.median(p3)
    p_vs_neg = sps.mannwhitneyu(p1, p0).pvalue
    return {
        "holds": bool(m0 <= m1 <= m3),
        "median_neg": float(m0), "median_int": float(m1), "median_pos": float(m3),
        "p_int_vs_neg": float(p_vs_neg),
        "n0": int(len(p0)), "n1": int(len(p1)), "n3": int(len(p3)),
    }


def _cross_slice_reproducibility(m, outcome):
    """② 跨切片复现性系数：在 CN/MCI/AD 各诊断层内检验「严格中间态」是否复现。

    返回 (repro_coeff, per_slice)：repro_coeff = 有效层中严格中间态成立的比例（[0,1]）。
    """
    layers = {}
    n_valid = 0
    n_holds = 0
    for dxv in [1.0, 2.0, 3.0]:
        layer = m[m.DX_bl == dxv]
        chk = _intermediate_check(layer, outcome)
        layers[DX_LABEL[dxv]] = chk
        if chk["holds"] is not None:
            n_valid += 1
            n_holds += int(chk["holds"])
    coeff = (n_holds / n_valid) if n_valid else None
    return {"coeff": coeff, "n_valid_slices": n_valid, "n_holds": n_holds, "per_slice": layers}


# --------------------------------------------------------------------------
# 探索环境主体
# --------------------------------------------------------------------------
class ExplorationEnvironment:
    def __init__(self, seed: int = 2026):
        self.seed = seed
        self.wide = load_wide()
        self.logs: list = []
        self.action_counter = 0
        self._aligned_cache: dict[str, pd.DataFrame] = {}

    # ---- 观察 ----
    def observe(self) -> dict:
        """返回 Agent 可读的当前状态（Page 2 §2.2 观察字段 + 环境元信息）。"""
        df = self.wide
        # 正确口径（exclude）：缺失血浆排除；四组计数与动作层一致
        plasma_ok = compute_plasma_status(df, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"], "exclude")
        correct = df.dropna(subset=["PET_STATUS"]).copy()
        correct["_PS"] = plasma_ok.loc[correct.index]
        correct = correct.dropna(subset=["_PS"])
        grp = correct["PET_STATUS"].astype(int) * 2 + correct["_PS"].astype(int)
        # 遗留列（subjects_wide.csv 的 PLASMA_STATUS，preprocess.py 曾把缺失编码为阴性）
        legacy = df.dropna(subset=["PET_STATUS", "PLASMA_STATUS"])
        lgrp = legacy["PET_STATUS"].astype(int) * 2 + legacy["PLASMA_STATUS"].astype(int)
        counts = {
            "n_total": int(len(df)),
            "n_with_marker_measured": int(len(correct)),
            "four_group_counts": {GROUP_LABELS[g]: int((grp == g).sum()) for g in GROUP_LABELS},
            "legacy_plasma_status_counts": {GROUP_LABELS[g]: int((lgrp == g).sum()) for g in GROUP_LABELS},
            "note": "four_group_counts=缺失血浆排除的正确口径；legacy_plasma_status_counts=初赛 PLASMA_STATUS 列（缺失编码为阴性）",
        }
        return {
            "observation": {
                "available_fields": ["AGE", "GENDER", "EDUCAT", "CDRSB_bl", "MMSE_bl", "ADAS13_bl",
                                     "DX_bl", "PTAU217_bl", "AB_RATIO_bl", "PET_STATUS", "PLASMA_STATUS",
                                     "CARDIO", "ENDO", "CKD", "D_ADAS13_yr", "D_CDRSB_yr"],
                "available_actions": list(ACTION_SCHEMAS.keys()),
                "action_schemas": ACTION_SCHEMAS,
            },
            "state": counts,
            "thresholds": {k: v for k, v in C.THRESHOLDS.items()},
            "rounds_done": len(self.logs),
        }

    def _validate_action_params(self, action: str, params: dict) -> None:
        if action == "define_discordance":
            if params.get("assay_pair", "ptau217_vs_pet") not in ACTION_SCHEMAS[action]["assay_pair"]:
                raise ValueError("非法 assay_pair")
            if params.get("outcome", "ADAS13") not in ACTION_SCHEMAS[action]["outcome"]:
                raise ValueError("非法 outcome")
            if params.get("missing_policy", "exclude") not in ACTION_SCHEMAS[action]["missing_policy"]:
                raise ValueError("非法 missing_policy")
        elif action == "discover_subtypes":
            if params.get("method", "kmeans") not in ACTION_SCHEMAS[action]["method"]:
                raise ValueError("非法 method")
            if int(params.get("k", 2)) not in ACTION_SCHEMAS[action]["k"]:
                raise ValueError("非法 k")
        elif action == "select_slice":
            if params.get("diagnosis", "CN") not in ACTION_SCHEMAS[action]["diagnosis"]:
                raise ValueError("非法 diagnosis")
            if params.get("age_band", "ALL") not in ACTION_SCHEMAS[action]["age_band"]:
                raise ValueError("非法 age_band")
        elif action == "sensitivity_analysis":
            if params.get("outcome", "CDRSB") not in ACTION_SCHEMAS[action]["outcome"]:
                raise ValueError("非法 outcome")
            if params.get("window", "all") not in ACTION_SCHEMAS[action]["window"]:
                raise ValueError("非法 window")
        elif action == "profile_mechanism":
            if params.get("target_group", "PET−/Plasma+") not in ACTION_SCHEMAS[action]["target_group"]:
                raise ValueError("非法 target_group")
            if params.get("alignment_window", "180d") not in ACTION_SCHEMAS[action]["alignment_window"]:
                raise ValueError("非法 alignment_window")
            if params.get("outcome", "ADAS13") not in ACTION_SCHEMAS[action]["outcome"]:
                raise ValueError("非法 outcome")
            markers = params.get("markers", ["tau_pet", "GFAP", "NfL"])
            if not isinstance(markers, list) or not set(markers).issubset({"tau_pet", "GFAP", "NfL"}):
                raise ValueError("非法 markers")

    # ---- 行动分发 ----
    def act(self, action: str, params: dict | None = None) -> dict:
        params = dict(params or {})
        self.action_counter += 1
        if action not in ACTION_SCHEMAS:
            return {"error": f"未知动作 {action}", "available": list(ACTION_SCHEMAS.keys())}
        try:
            self._validate_action_params(action, params)
            fn = getattr(self, f"_action_{action}")
            fb = fn(params)
        except Exception as e:  # 失败轮次也返回结构化反馈，不中断探索闭环
            fb = {"error": f"{type(e).__name__}: {e}"}
        fb["action"] = action
        fb["params"] = params
        return _clean(fb)

    # ---- 记录 ----
    def record(self, round_id: int, action: str, params: dict, rationale: str,
               feedback: dict, state_update: str, success: bool | None = None) -> dict:
        entry = {
            "round_id": int(round_id),
            "timestamp": _now(),
            "action": action,
            "params": params,
            "rationale": rationale,
            "state_update": state_update,  # 定义 / 切片 / 方法 / 假设 之一（Page 3 §3.3 更新约定）
            "success": success,
            "feedback": _clean(feedback),
        }
        self.logs.append(entry)
        return entry

    def flush_logs(self) -> Path:
        """把当前 logs 写回 data/processed/exploration_log.jsonl（幂等，覆盖写）。"""
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_PATH, "w", encoding="utf-8") as f:
            for e in self.logs:
                f.write(json.dumps(e, ensure_ascii=False) + "\n")
        return LOG_PATH

    # ---- 动作 1：define_discordance ----
    def _action_define_discordance(self, p: dict) -> dict:
        assay = p.get("assay_pair", "ptau217_vs_pet")
        outcome = OUTCOME_COL[p.get("outcome", "ADAS13")]
        missing_policy = p.get("missing_policy", "exclude")
        if assay == "ptau217_vs_pet":
            thr = float(p.get("plasma_threshold", C.THRESHOLDS["p_tau217_pg_ml"]))
            marker_col = "PTAU217_bl"
        else:
            thr = float(p.get("plasma_threshold", C.THRESHOLDS["abeta_ratio"]))
            marker_col = "AB_RATIO_bl"

        plasma = compute_plasma_status(self.wide, assay, thr, missing_policy)
        sub = add_group(self.wide, plasma)
        m = build_longitudinal(sub)
        m = make_comorbidity(m)

        # ① 不一致组规模与构成
        n_by_group = {GROUP_LABELS[g]: int((m.GROUP == g).sum()) for g in GROUP_LABELS}
        n_all = int(len(m))
        n_discord = n_by_group[GROUP_LABELS[1]] + n_by_group[GROUP_LABELS[2]]
        # 缺失标志物按 PET 状态分层（暴露「未测血浆被当成阴性」的膨胀量）
        pet_only = self.wide.dropna(subset=["PET_STATUS"])
        marker_missing_by_pet = {int(k): int(v) for k, v in
                                 pet_only.groupby("PET_STATUS")[marker_col].apply(lambda s: int(s.isna().sum())).items()}
        size = {"n_total": n_all, "n_discordant": n_discord,
                "discordance_rate": (n_discord / n_all if n_all else None),
                "by_group": n_by_group,
                "missing_policy": missing_policy,
                "marker_missing_by_pet_status": marker_missing_by_pet,
                "composition": {"age": {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, "AGE"].median()) for g in GROUP_LABELS},
                                "female_frac": {GROUP_LABELS[g]: float((m.loc[m.GROUP == g, "GENDER"] == 2).mean()) for g in GROUP_LABELS}}}

        # ② 跨切片复现性
        repro = _cross_slice_reproducibility(m, outcome)

        # ③ 轨迹差异（中间态 vs 双阴 / 双阳）
        p_neg, n1, n0 = _mannwhitney(m, "GROUP", outcome, 1, 0)
        p_pos, _, n3 = _mannwhitney(m, "GROUP", outcome, 1, 3)
        es_neg = _cohens_d_ci(m.loc[m.GROUP == 1, outcome], m.loc[m.GROUP == 0, outcome])
        med = {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, outcome].median()) for g in GROUP_LABELS}
        traj = {"median_by_group": med, "intermediate_vs_neg": {"p": p_neg, "cohens_d": es_neg},
                "intermediate_vs_pos": {"p": p_pos}}

        # ④ 混杂控制后效应变化
        base = _adjusted_coeffs(m, "GROUP", outcome, ["AGE", "GENDER", "EDUCAT"], ref_val=0)
        full = _adjusted_coeffs(m, "GROUP", outcome, ["AGE", "GENDER", "EDUCAT", "CARDIO", "ENDO", "CKD"], ref_val=0)
        conf = {"base_coefs": base.get("group_coefs") if base else None,
                "comorbidity_adjusted_coefs": full.get("group_coefs") if full else None}

        # ⑤ 与随机参照差异
        H_obs, p_perm = _permutation_p(m, "GROUP", outcome, n_perm=500, seed=self.seed)

        return {"metrics": {
            "discordant_group_size": size,
            "cross_slice_reproducibility": repro,
            "trajectory_separation": traj,
            "confounder_adjusted_change": conf,
            "vs_random_reference": {"H_obs": H_obs, "perm_p": p_perm},
        }, "interpretation": _interpret_discordance(repro, traj, p_perm)}

    # ---- 动作 2：discover_subtypes ----
    def _action_discover_subtypes(self, p: dict) -> dict:
        method = p.get("method", "kmeans")
        k = int(p.get("k", 2))
        diagnosis = p.get("diagnosis", "ALL")
        features = p.get("features", ["PTAU217_bl", "AB_RATIO_bl", "CARDIO", "ENDO"])
        outcome = OUTCOME_COL["ADAS13"]

        plasma = compute_plasma_status(self.wide, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"])
        sub = add_group(self.wide, plasma)
        m = build_longitudinal(sub)
        m = make_comorbidity(m)
        pm = m[(m.PET_STATUS == 0) & (m.PLASMA_STATUS == 1)].copy()
        if diagnosis != "ALL":
            pm = pm[pm.DX_bl == {"CN": 1.0, "MCI": 2.0, "AD": 3.0}[diagnosis]]
        pm = pm.dropna(subset=[f for f in features if f in pm.columns])

        X = StandardScaler().fit_transform(pm[features])
        if method == "kmeans":
            model = KMeans(n_clusters=k, n_init=10, random_state=42).fit(X)
            labels = model.labels_
            sil = silhouette_score(X, labels)
        elif method == "gmm":
            model = GaussianMixture(n_components=k, random_state=42).fit(X)
            labels = model.predict(X)
            sil = silhouette_score(X, labels)
        else:
            raise ValueError(f"未知 method: {method}")
        pm["CLUSTER"] = labels

        # ① 规模与构成
        size = {"diagnosis": diagnosis, "n_target_group": int(len(pm)),
                "cluster_sizes": {int(c): int((labels == c).sum()) for c in np.unique(labels)},
                "silhouette": float(sil),
                "stability": _cluster_stability(X, labels, method, k, seed=self.seed)}

        # ② 跨切片复现（CN vs MCI 层内重聚类，看合并症主导的簇是否稳定出现）
        repro = {}
        for dxv, dname in [(1.0, "CN"), (2.0, "MCI")]:
            layer = pm[pm.DX_bl == dxv]
            if len(layer) < 30:
                repro[dname] = {"reproduced": None, "n": int(len(layer))}
                continue
            Xl = StandardScaler().fit_transform(layer[features])
            lab = KMeans(n_clusters=k, n_init=10, random_state=42).fit_predict(Xl) if method == "kmeans" else \
                  GaussianMixture(n_components=k, random_state=42).fit(Xl).predict(Xl)
            # 简化：合并症负担最高的簇是否稳定存在（双簇 comorbidity 差）
            layer = layer.copy(); layer["CL"] = lab
            com_diff = layer.groupby("CL")["CARDIO"].mean().max() - layer.groupby("CL")["CARDIO"].mean().min()
            repro[dname] = {"reproduced": bool(com_diff > 0.2), "n": int(len(layer)),
                            "comorbidity_split": float(com_diff)}

        # ③ 轨迹分离（簇间）
        cl = sorted(pm.CLUSTER.unique())
        traj = {}
        if len(cl) >= 2:
            p_mw = sps.mannwhitneyu(pm.loc[pm.CLUSTER == cl[0], outcome].dropna(),
                                    pm.loc[pm.CLUSTER == cl[1], outcome].dropna()).pvalue
            es = _cohens_d_ci(pm.loc[pm.CLUSTER == cl[0], outcome], pm.loc[pm.CLUSTER == cl[1], outcome])
            traj = {"cluster_medians": {int(c): float(pm.loc[pm.CLUSTER == c, outcome].median()) for c in cl},
                    "mw_p": float(p_mw), "cohens_d": es}

        # ④ 混杂控制后效应变化（簇 + AGE/GENDER/EDUCAT）
        pm = pm.dropna(subset=["AGE", "GENDER", "EDUCAT", outcome])
        if len(pm) >= 20 and len(pm.CLUSTER.unique()) >= 2:
            Xr = pd.DataFrame({"C": (pm.CLUSTER == cl[0]).astype(int), "AGE": pm.AGE,
                               "GENDER": pm.GENDER.map({1: 1, 2: 0}), "EDUCAT": pm.EDUCAT})
            yr = pm.loc[Xr.index, outcome]
            lm = LinearRegression().fit(Xr, yr)
            conf = {"cluster_coef": float(lm.coef_[0]), "cluster_coef_adjusted_ci_note": "系数=簇差（Δ/年），校正年龄/性别/教育后"}
        else:
            conf = {"cluster_coef": None}

        # ⑤ 与随机参照差异（打乱簇标签重算轨迹分离）
        H_obs, p_perm = _permutation_p(pm, "CLUSTER", outcome, n_perm=300, seed=self.seed)

        return {"metrics": {
            "discordant_group_size": size,
            "cross_slice_reproducibility": repro,
            "trajectory_separation": traj,
            "confounder_adjusted_change": conf,
            "vs_random_reference": {"H_obs": H_obs, "perm_p": p_perm},
        }, "interpretation": _interpret_subtypes(sil, traj, p_perm, size)}

    # ---- 动作 3：select_slice ----
    def _action_select_slice(self, p: dict) -> dict:
        diag = p.get("diagnosis", "CN")
        age_band = p.get("age_band", "ALL")
        outcome = OUTCOME_COL["ADAS13"]

        plasma = compute_plasma_status(self.wide, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"])
        sub = add_group(self.wide, plasma)
        m = build_longitudinal(sub)
        if diag != "ALL":
            m = m[m.DX_bl == {"CN": 1.0, "MCI": 2.0, "AD": 3.0}[diag]]
        if age_band == "<75":
            m = m[m.AGE < 75]
        elif age_band == ">=75":
            m = m[m.AGE >= 75]

        size = {"n_slice": int(len(m)),
                "by_group": {GROUP_LABELS[g]: int((m.GROUP == g).sum()) for g in GROUP_LABELS}}
        repro = _intermediate_check(m, outcome)
        med = {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, outcome].median()) for g in GROUP_LABELS}
        p_neg, _, _ = _mannwhitney(m, "GROUP", outcome, 1, 0)
        p_pos, _, _ = _mannwhitney(m, "GROUP", outcome, 1, 3)
        es_neg = _cohens_d_ci(m.loc[m.GROUP == 1, outcome], m.loc[m.GROUP == 0, outcome])
        conf = _adjusted_coeffs(m, "GROUP", outcome, ["AGE", "GENDER", "EDUCAT"], ref_val=0)
        H_obs, p_perm = _permutation_p(m, "GROUP", outcome, n_perm=300, seed=self.seed)

        return {"metrics": {
            "discordant_group_size": size,
            "cross_slice_reproducibility": {"slice_intermediate_check": repro,
                                             "note": "本动作在指定切片内检验中间态是否成立"},
            "trajectory_separation": {"median_by_group": med,
                                       "intermediate_vs_neg": {"p": p_neg, "cohens_d": es_neg},
                                       "intermediate_vs_pos": {"p": p_pos}},
            "confounder_adjusted_change": conf.get("group_coefs") if conf else None,
            "vs_random_reference": {"H_obs": H_obs, "perm_p": p_perm},
        }, "interpretation": _interpret_slice(diag, age_band, repro, med)}

    # ---- 动作 4：test_confounder ----
    def _action_test_confounder(self, p: dict) -> dict:
        control_vars = p.get("control_vars", ["CARDIO", "ENDO"])
        outcome = OUTCOME_COL["ADAS13"]

        plasma = compute_plasma_status(self.wide, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"])
        sub = add_group(self.wide, plasma)
        m = build_longitudinal(sub)
        m = make_comorbidity(m)
        m["GENDER"] = m["GENDER"].map({1: 1, 2: 0})

        # 逐级加入控制变量，看中间态（GROUP==1）系数如何变化
        base_vars = ["AGE", "GENDER", "EDUCAT"]
        size = {"n_target": int((m.GROUP == 1).sum())}
        ladder = {}
        for extra in [[], control_vars]:
            varset = base_vars + extra
            d = m.dropna(subset=["GROUP", outcome] + varset)
            if d.empty:
                ladder["+".join(varset)] = None
                continue
            gcols = {f"G{g}": (d.GROUP.astype(int) == g).astype(int) for g in [1, 2, 3]}
            X = pd.DataFrame(gcols)
            for cv in varset:
                X[cv] = pd.to_numeric(d[cv], errors="coerce")
            X = X.dropna()
            y = d.loc[X.index, outcome]
            lm = LinearRegression().fit(X, y)
            key = "base" if not extra else "+".join(extra)
            ladder[key] = {"coef_G1": float(lm.coef_[0]), "coef_G2": float(lm.coef_[1]),
                           "coef_G3": float(lm.coef_[2]), "r2": float(lm.score(X, y))}

        base_coef = ladder.get("base", {}).get("coef_G1")
        full_coef = ladder.get("+".join(control_vars), {}).get("coef_G1") if control_vars else base_coef
        delta = (full_coef - base_coef) if (base_coef is not None and full_coef is not None) else None

        # ③④⑤：以最终模型为基准做轨迹分离与随机参照
        med = {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, outcome].median()) for g in GROUP_LABELS}
        p_neg, _, _ = _mannwhitney(m, "GROUP", outcome, 1, 0)
        es_neg = _cohens_d_ci(m.loc[m.GROUP == 1, outcome], m.loc[m.GROUP == 0, outcome])
        H_obs, p_perm = _permutation_p(m, "GROUP", outcome, n_perm=300, seed=self.seed)

        return {"metrics": {
            "discordant_group_size": size,
            "cross_slice_reproducibility": {"confounder_ladder": ladder,
                                             "note": "中间态系数随控制变量加入的稳定性（定义敏感性）"},
            "trajectory_separation": {"median_by_group": med,
                                       "intermediate_vs_neg": {"p": p_neg, "cohens_d": es_neg}},
            "confounder_adjusted_change": {"base_coef_G1": base_coef, "adjusted_coef_G1": full_coef,
                                            "delta": delta, "control_vars": control_vars},
            "vs_random_reference": {"H_obs": H_obs, "perm_p": p_perm},
        }, "interpretation": _interpret_confounder(delta, p_neg)}

    # ---- 动作 5：sensitivity_analysis ----
    def _action_sensitivity_analysis(self, p: dict) -> dict:
        outcome_key = p.get("outcome", "CDRSB")
        window = p.get("window", "all")
        outcome = OUTCOME_COL[outcome_key]

        plasma = compute_plasma_status(self.wide, "ptau217_vs_pet", C.THRESHOLDS["p_tau217_pg_ml"])
        sub = add_group(self.wide, plasma)
        m = build_longitudinal(sub)
        if window == ">=2yr_followup":
            yrs_col = "YRS_ADAS" if outcome_key == "ADAS13" else "YRS_CDR"
            m = m[m[yrs_col] >= 2]

        size = {"n_total": int(len(m)),
                "by_group": {GROUP_LABELS[g]: int((m.GROUP == g).sum()) for g in GROUP_LABELS}}
        repro = _cross_slice_reproducibility(m, outcome)
        med = {GROUP_LABELS[g]: float(m.loc[m.GROUP == g, outcome].median()) for g in GROUP_LABELS}
        p_neg, _, _ = _mannwhitney(m, "GROUP", outcome, 1, 0)
        p_pos, _, _ = _mannwhitney(m, "GROUP", outcome, 1, 3)
        es_neg = _cohens_d_ci(m.loc[m.GROUP == 1, outcome], m.loc[m.GROUP == 0, outcome])
        conf = _adjusted_coeffs(m, "GROUP", outcome, ["AGE", "GENDER", "EDUCAT"], ref_val=0)
        H_obs, p_perm = _permutation_p(m, "GROUP", outcome, n_perm=300, seed=self.seed)

        return {"metrics": {
            "discordant_group_size": size,
            "cross_slice_reproducibility": repro,
            "trajectory_separation": {"median_by_group": med,
                                       "intermediate_vs_neg": {"p": p_neg, "cohens_d": es_neg},
                                       "intermediate_vs_pos": {"p": p_pos}},
            "confounder_adjusted_change": conf.get("group_coefs") if conf else None,
            "vs_random_reference": {"H_obs": H_obs, "perm_p": p_perm},
        }, "interpretation": _interpret_sensitivity(outcome_key, window, repro, med, p_neg)}


    # ---- 动作 6：profile_mechanism ----
    def _action_profile_mechanism(self, p: dict) -> dict:
        from mechanism_analysis import profile_mechanism
        target_group = p.get("target_group", "PET−/Plasma+")
        markers = p.get("markers", ["tau_pet", "GFAP", "NfL"])
        window = p.get("alignment_window", "180d")
        outcome = p.get("outcome", "ADAS13")
        if window not in ALIGNMENT_WINDOWS:
            raise ValueError(f"alignment_window 必须是 {sorted(ALIGNMENT_WINDOWS)}")
        if window not in self._aligned_cache:
            from data_layer import build_aligned_cohort
            self._aligned_cache[window] = build_aligned_cohort(window=window, missing_policy="exclude")
        frame = self._aligned_cache[window]
        profile = profile_mechanism(frame, target_group=target_group, markers=markers,
                                    outcome=outcome, alignment_window=window, seed=self.seed)
        ranking = profile.get("candidate_mechanism_ranking", [])
        return {
            "metrics": {
                "discordant_group_size": {"n_aligned": int(len(frame)),
                    "by_group": {GROUP_LABELS[g]: int((frame.GROUP == g).sum()) for g in GROUP_LABELS}},
                "cross_slice_reproducibility": {"alignment_window": window,
                    "time_gap_days": profile["lag_evidence"].get("time_gap_days")},
                "trajectory_separation": profile["trajectory_evidence"],
                "confounder_adjusted_change": profile["trajectory_evidence"].get("adjusted"),
                "vs_random_reference": {"multiple_testing": profile.get("multiple_testing")},
                "mechanism_evidence": {
                    "lag": profile["lag_evidence"], "injury": profile["injury_evidence"],
                    "tau": profile["tau_evidence"], "trajectory": profile["trajectory_evidence"],
                    "missingness": profile["missingness_evidence"],
                },
                "candidate_mechanisms": ranking,
                "negative_evidence": profile.get("negative_evidence", []),
                "next_best_action": profile.get("next_best_action"),
            },
            "interpretation": "；".join(f"{item['mechanism']}={item['score']:.2f}" for item in ranking),
        }


# --------------------------------------------------------------------------
# 解读辅助（Agent 的"受挫/修正"依据，也写进日志 rationale 由 Agent 生成）
# --------------------------------------------------------------------------
def _interpret_discordance(repro, traj, p_perm):
    holds = repro.get("coeff")
    p1 = traj["intermediate_vs_neg"].get("p")
    s = []
    s.append(f"不一致率见规模指标；中间态轨迹 vs 双阴 p={p1 if p1 is None else round(p1,4)}")
    s.append(f"跨切片复现系数={holds if holds is None else round(holds,2)}（CN/MCI/AD 层内严格中间态比例）")
    s.append(f"置换检验 p={None if p_perm is None else round(p_perm,4)}（<0.05 视为超随机）")
    return "；".join(s)


def _interpret_subtypes(sil, traj, p_perm, size):
    p_mw = traj.get("mw_p")
    s = [f"轮廓系数={round(sil,3)}"]
    if p_mw is not None:
        s.append(f"两簇轨迹 p={round(p_mw,4)}（≥0.05 → 轨迹不可分，负结果）")
    s.append(f"置换 p={None if p_perm is None else round(p_perm,4)}")
    return "；".join(s)


def _interpret_slice(diag, age_band, repro, med):
    holds = repro.get("holds")
    return f"{diag}/{age_band} 层：严格中间态={'成立' if holds else ('不成立' if holds is False else '样本不足')}；" \
           f"双阴中位={med.get('PET−/Plasma−')} 中间={med.get('PET−/Plasma+')} 双阳={med.get('PET+/Plasma+')}"


def _interpret_confounder(delta, p_neg):
    d = None if delta is None else round(delta, 4)
    return f"加入合并症后中间态系数变化={d}；中间态 vs 双阴 p={None if p_neg is None else round(p_neg,4)}"


def _interpret_sensitivity(outcome_key, window, repro, med, p_neg):
    return f"{outcome_key}/{window}：复现系数={None if repro.get('coeff') is None else round(repro['coeff'],2)}；" \
           f"中间态 vs 双阴 p={None if p_neg is None else round(p_neg,4)}"


if __name__ == "__main__":
    # 冒烟测试：跑一次 observe + 一个动作
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    env = ExplorationEnvironment()
    obs = env.observe()
    print(json.dumps(obs, ensure_ascii=False, indent=1)[:2000])
    fb = env.act("define_discordance", {"assay_pair": "ptau217_vs_pet"})
    print(json.dumps(fb, ensure_ascii=False, indent=1)[:2000])
