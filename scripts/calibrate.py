# -*- coding: utf-8 -*-
"""calibrate.py — 以 Aβ PET 官方状态为参考，标定血浆标志物阳性阈值（Youden J）

用法：python scripts/calibrate.py
输出：data/processed/calibrate_report.txt（AUC、Youden 阈值、四格表、建议值）

⚠️ 方法透明性：以 PET 为参考标定血浆阈值，与 E1/Pyun 2024 的做法同类；
该阈值只定义"不一致"的切分，不构成独立诊断结论。阈值由人确认后填入 config。
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import config as C


def youden(score, y):
    """Youden J 最优阈值。返回 (threshold, J, spec, sens)。"""
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y, score)
    j = tpr - fpr
    i = int(np.argmax(j))
    return float(thr[i]), float(j[i]), float(1 - fpr[i]), float(tpr[i])


def main():
    df = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    sub = df.dropna(subset=["PTAU217_bl", "PET_STATUS"]).copy()
    lines = [f"校准样本：{len(sub)} 人（血浆 pT217_F × PET 状态均有）"]
    lines.append(f"PET 阳性率：{sub.PET_STATUS.mean():.3f}")

    report = {}
    # direction: +1 正向标志物（越高越病理），-1 逆向（越低越病理，如 Aβ42/40 比值）
    for col, key, name, direction in [
        ("PTAU217_bl", "p_tau217_pg_ml", "pT217_F", +1),
        ("AB_RATIO_bl", "abeta_ratio", "AB42/40", -1),
    ]:
        # 每个标志物单独 dropna：AB42/40 与 pT217 的缺失模式不完全重合
        sub2 = sub.dropna(subset=[col])
        s = sub2[col].astype(float)
        y = sub2["PET_STATUS"].astype(int)
        from sklearn.metrics import roc_auc_score
        auc = roc_auc_score(y, s * direction)
        thr, j, spec, sens = youden(s * direction, y)
        thr = thr / direction  # 还原到原始刻度
        # 四格表（按该阈值；逆向标志物为 s < thr）
        pred = ((s > thr) if direction == +1 else (s < thr)).astype(int)
        tp = ((pred == 1) & (y == 1)).sum(); fp = ((pred == 1) & (y == 0)).sum()
        tn = ((pred == 0) & (y == 0)).sum(); fn = ((pred == 0) & (y == 1)).sum()
        n_discord = fp + fn
        lines += [
            f"\n===== {name} =====",
            f"AUC = {auc:.3f}",
            f"Youden 阈值 = {thr:.3f}  (J={j:.3f}, sens={sens:.2f}, spec={spec:.2f})",
            f"四格表: TP={tp} FP={fp} FN={fn} TN={tn}  → 不一致 {n_discord} 例 ({n_discord/len(sub2):.1%})",
        ]
        report[key] = {"auc": round(auc, 3), "youden_thr": round(thr, 3), "n_discord": int(n_discord)}

    # 建议阈值（保守取整，减少定义敏感性）
    lines.append("\n===== 建议写入 config.THRESHOLDS 的值 =====")
    for key in ("p_tau217_pg_ml", "abeta_ratio"):
        if report.get(key):
            lines.append(f"{key} = {report[key]['youden_thr']}")
    txt = "\n".join(lines)
    print(txt)
    (C.PROC_DIR / "calibrate_report.txt").write_text(txt, encoding="utf-8")
    (C.PROC_DIR / "calibrate.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
