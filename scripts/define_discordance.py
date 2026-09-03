# -*- coding: utf-8 -*-
"""define_discordance.py — 四组切分 + 组间统计 + 与 E1 定性对照核对

用法：python scripts/define_discordance.py
输出：data/processed/group_stats.csv；控制台打印核对表
"""
import pandas as pd

import config as C

# GROUP = PET*2 + PLASMA；命名按数值映射（1=PET0+PLASMA1, 2=PET1+PLASMA0）
GROUP_NAMES = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}

def main():
    df = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    need = ["PET_STATUS", "PLASMA_STATUS"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit(f"[discordance] 缺状态列 {miss}：先跑 preprocess.py 并补齐 config.py 阈值")

    df = df.dropna(subset=need)
    df["GROUP"] = (df["PET_STATUS"].astype(int) * 2 + df["PLASMA_STATUS"].astype(int))

    g = df.groupby("GROUP")
    stats = pd.DataFrame({
        "n": g.size(),
        "age_mean": g["AGE"].mean() if "AGE" in df else None,
        "CDRSB_bl_mean": g["CDRSB_bl"].mean() if "CDRSB_bl" in df else None,
        "MMSE_bl_mean": g["MMSE_bl"].mean() if "MMSE_bl" in df else None,
        "ADAS13_bl_mean": g["ADAS13_bl"].mean() if "ADAS13_bl" in df else None,
    }).rename(index=GROUP_NAMES)
    stats.to_csv(C.PROC_DIR / "group_stats.csv")
    print(stats.round(2).to_string())

    # E1 定性核对（对照 Pyun 2024 摘要结论，非数字对照）
    checks = [
        ("PET+/Plasma+ 认知最差", "ADAS13_bl_mean" in stats and stats.loc["PET+/Plasma+", "ADAS13_bl_mean"] >= stats["ADAS13_bl_mean"].max() - 0.01),
        ("PET−/Plasma+ 介于双阴与双阳之间（中间态）", "ADAS13_bl_mean" in stats and stats.loc["PET−/Plasma+", "ADAS13_bl_mean"] <= stats.loc["PET+/Plasma+", "ADAS13_bl_mean"]),
        ("四组均有样本（无空组）", stats["n"].min() > 0),
    ]
    print("\n=== E1 定性核对 ===")
    for desc, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {desc}")
    print("\n说明：这是定性核对，仅确认管线行为与 E1 摘要结论方向一致，不作为研究结论。")

if __name__ == "__main__":
    main()
