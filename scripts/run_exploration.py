# -*- coding: utf-8 -*-
"""run_exploration.py — 一键复现复赛三件套（探索环境 + 参照系 + 探索日志）

用法（在 goai-open-exploration/ 下）：
    python scripts/run_exploration.py              # MockPolicy（无 key，确定性）
    LLM_API_KEY=... python scripts/run_exploration.py   # LLM 驱动闭环

输出（全部落盘 data/processed/）：
    exploration_log.jsonl   逐轮探索日志（round_log schema，Page 2 §2.3）
    exploration_report.md   人可读探索报告（评审友好）
    reference_report.txt/.json  四层参照系（Page 3 §3.2）
"""
import json
import sys
from pathlib import Path

import config as C
from environment import ExplorationEnvironment
from reference import ReferenceFrame
from agent import ExplorationAgent
from preprocess import main as preprocess_main
from mechanism_analysis import run_mechanism_analysis
from data_layer import validate_raw_inputs
from make_dashboard import main as make_dashboard_main
from validation import run_validation


def write_markdown_report(logs: list, path: Path) -> Path:
    lines = ["# 探索报告（GOAI AI4S 探索赛道 · 复赛）\n",
             "> 由 `scripts/run_exploration.py` 自动生成；逐轮可复现，日志见 `exploration_log.jsonl`。\n"]
    lines.append("| 轮 | 动作 | 参数 | 状态更新 | 关键结果 |")
    lines.append("|---|---|---|---|---|")
    def _f(v, nd=3):
        return "—" if v is None else f"{v:.{nd}g}"
    for e in logs:
        fb = e["feedback"]
        m = fb.get("metrics", {})
        key = ""
        if e["action"] == "define_discordance":
            s = m["discordant_group_size"]
            key = f"不一致率 {s['discordance_rate']:.1%}；四组 {s['by_group']}"
        elif e["action"] == "select_slice":
            chk = m["cross_slice_reproducibility"]["slice_intermediate_check"]
            key = f"中间态={'成立' if chk['holds'] else '不成立'}（p_int_vs_neg={_f(chk.get('p_int_vs_neg'))}）"
        elif e["action"] == "test_confounder":
            key = f"Δ系数={_f(m['confounder_adjusted_change'].get('delta'), 4)}"
        elif e["action"] == "discover_subtypes":
            key = f"轮廓={m['discordant_group_size']['silhouette']:.2f}；轨迹 p={_f(m['trajectory_separation'].get('mw_p'))}"
        elif e["action"] == "sensitivity_analysis":
            key = f"复现系数={_f(m['cross_slice_reproducibility'].get('coeff'))}"
        elif e["action"] == "profile_mechanism":
            ranking = m.get("candidate_mechanisms", [])
            key = "机制排序=" + ", ".join(item.get("mechanism", "?") for item in ranking[:3])
        lines.append(f"| {e['round_id']} | {e['action']} | `{json.dumps(e['params'], ensure_ascii=False)}` "
                     f"| {e['state_update']} | {key} |")
    lines.append("")
    lines.append("## 每轮 rationale")
    lines.append("")
    for e in logs:
        lines.append(f"- **Round {e['round_id']}**（{e['action']}，{e['state_update']}）：{e['rationale']}")
    lines.append("")
    txt = "\n".join(lines)
    path.write_text(txt, encoding="utf-8")
    return path


def main():
    n_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    raw_report = validate_raw_inputs()
    (C.PROC_DIR / "raw_validation.json").write_text(json.dumps(raw_report, ensure_ascii=False, indent=2), encoding="utf-8")
    if not raw_report["valid"]:
        raise SystemExit("[pipeline] raw data schema validation failed")
    preprocess_main()
    mechanism_report = run_mechanism_analysis("180d")
    validation_report = run_validation("180d")
    env = ExplorationEnvironment()

    # ① 最小可运行探索环境：冒烟（observe + 全部 6 个动作各跑一次）
    obs = env.observe()
    print(f"[env] 观察：样本 {obs['state']['n_total']}，标志物已测 {obs['state']['n_with_marker_measured']}，"
          f"四组(修正口径) {obs['state']['four_group_counts']}")

    # ② 参照系
    ref = ReferenceFrame(env)
    ref_path = ref.write_report()
    print(f"[reference] 参照系报告已写：{ref_path}")

    # ③ Agent 探索闭环
    agent = ExplorationAgent(env)
    logs, log_path = agent.run(n_rounds=n_rounds)
    md_path = write_markdown_report(logs, C.PROC_DIR / "exploration_report.md")
    print(f"[agent] 探索日志已写：{log_path}（{len(logs)} 轮）")
    print(f"[agent] 探索报告已写：{md_path}")
    print(f"[mechanism] 对齐样本 {mechanism_report['n_aligned']}；报告已写：{C.PROC_DIR / 'mechanism_report.json'}")
    print(f"[validation] 验证样本 {validation_report['split']['n_validation']}；报告已写：{C.PROC_DIR / 'validation_report.json'}")
    make_dashboard_main()

    print("\n复赛三件套产出：")
    print("  ① 最小可运行探索环境  scripts/environment.py（observe/act/record + 6 动作 + 机制反馈）")
    print("  ② 参照系设计          scripts/reference.py → data/processed/reference_report.txt")
    print("  ③ 探索日志            data/processed/exploration_log.jsonl + exploration_report.md")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
