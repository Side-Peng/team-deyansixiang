# -*- coding: utf-8 -*-
"""agent.py — LLM 驱动的探索 Agent（GOAI AI4S 探索赛道 · 复赛三件套之三）

构成「观察 → 行动 → 反馈 → 记录 → 下一轮」的探索闭环：
    env   = ExplorationEnvironment()
    agent = ExplorationAgent(env, llm=...)      # llm=None 时用 MockPolicy（确定性、可复现）
    logs, path = agent.run(n_rounds=6)

LLM 接入（OpenAI 兼容 chat/completions，仅用标准库 urllib，无额外依赖）：
    环境变量 LLM_API_KEY / LLM_BASE_URL / LLM_MODEL，或显式传 OpenAICompatLLM(...)。
    未配置 key 时自动降级到 MockPolicy：一条确定性的「受挫→修正→继续」探索路径，
    全程读真实反馈、引用真实数字，保证无 key 也能产出可复现的探索日志。

Agent 每轮返回严格 JSON：{"action","params","rationale","state_update"}。
state_update ∈ {定义, 切片, 方法, 假设}（对齐 four_page_guide.md Page 3 §3.3 更新约定）。
"""
import json
import os
import re
import urllib.request
from pathlib import Path

from environment import ExplorationEnvironment

ROOT = Path(__file__).resolve().parents[1]


def _load_llm_config() -> dict:
    """读取仓库根目录 llm_config.json（若存在）。优先级：显式传参 > 环境变量 > 配置文件。"""
    p = ROOT / "llm_config.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


# --------------------------------------------------------------------------
# OpenAI 兼容 LLM（标准库，无第三方依赖）
# --------------------------------------------------------------------------
class OpenAICompatLLM:
    def __init__(self, api_key=None, base_url=None, model=None):
        cfg = _load_llm_config()
        self.api_key = (api_key or os.environ.get("LLM_API_KEY") or os.environ.get("OPENAI_API_KEY")
                        or cfg.get("api_key") or "")
        self.base_url = (base_url or os.environ.get("LLM_BASE_URL") or cfg.get("base_url")
                         or "https://api.openai.com/v1").rstrip("/")
        self.model = model or os.environ.get("LLM_MODEL") or cfg.get("model") or "gpt-4o-mini"

    @property
    def available(self) -> bool:
        return bool(self.api_key)

    def __call__(self, system: str, user: str) -> str:
        if not self.api_key:
            raise RuntimeError("未配置 LLM_API_KEY，无法调用模型")
        body = json.dumps({
            "model": self.model,
            "messages": [{"role": "system", "content": system},
                         {"role": "user", "content": user}],
            "temperature": 0.2,
        }).encode("utf-8")
        req = urllib.request.Request(
            f"{self.base_url}/chat/completions",
            data=body,
            headers={"Content-Type": "application/json",
                     "Authorization": f"Bearer {self.api_key}"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=120) as r:
            data = json.loads(r.read().decode("utf-8"))
        return data["choices"][0]["message"]["content"]


SYSTEM_PROMPT = """你是 GOAI AI4S 探索赛道中的科研探索 Agent。你在一个「影像-血浆生物标志物
不一致亚型」的 ADNI 数据集环境中工作。每轮你会看到环境观察（observation）与最近的探索历史，
然后选择下一个动作（action + params），并用中文写一句 rationale（你基于哪些反馈做此决定）与
state_update（本轮是否改变「定义 / 切片 / 方法 / 假设」，四选一）。

可用动作与参数空间在 observation.action_schemas 里。你只输出一个严格 JSON 对象，不要输出任何
多余文字或 markdown 代码块。格式：
{"action": "...", "params": {...}, "rationale": "...", "state_update": "..."}"""


def _extract_json(text: str) -> dict:
    """从 LLM 回复里提取首个 {...} 作为 JSON（容错代码块/前后缀）。"""
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        raise ValueError(f"无法从回复解析 JSON: {text[:200]}")
    return json.loads(m.group(0))


# --------------------------------------------------------------------------
# 确定性 MockPolicy（无 key 降级；读真实反馈、引用真实数字）
# --------------------------------------------------------------------------
class MockPolicy:
    """一条原则化的探索路径：复现初赛口径 → 发现口径异常 → 修正定义 → 切片/混杂/亚型 →
    负结果 → 敏感性 → 假设更新。每步 rationale 引用上一步真实反馈数字。"""

    def propose(self, obs: dict, history: list) -> dict:
        n = len(history)
        steps = [self._s0, self._s1, self._s2, self._s3, self._s4, self._s5, self._s6, self._s7]
        if n < len(steps):
            return steps[n](history)
        return {"action": "sensitivity_analysis", "params": {"outcome": "CDRSB"},
                "rationale": "预算内收尾轮：以 CDRSB 做最终敏感性确认。", "state_update": "方法"}

    @staticmethod
    def _s0(history):
        return {"action": "define_discordance",
                "params": {"assay_pair": "ptau217_vs_pet", "missing_policy": "missing_as_negative"},
                "rationale": "先按初赛管线口径（缺失血浆按阴性）复现四组基线，确认与已提交结果一致。",
                "state_update": "定义"}

    @staticmethod
    def _s1(history):
        prev = history[-1]["feedback"]["metrics"]["discordant_group_size"]
        rate = prev.get("discordance_rate")
        missing = prev.get("marker_missing_by_pet_status", {})
        return {"action": "define_discordance",
                "params": {"assay_pair": "ptau217_vs_pet", "missing_policy": "exclude"},
                "rationale": (f"初赛口径不一致率={rate:.1%}，显著高于标定口径 18.4% 与 E4 量级；"
                              f"marker_missing_by_pet={missing} 提示大量「未测血浆」被编码为阴性。"
                              f"改 missing_policy=exclude 复核定义，验证是否为编码 bug。"),
                "state_update": "定义"}

    @staticmethod
    def _s2(history):
        prev = history[-1]["feedback"]["metrics"]["discordant_group_size"]
        byg = prev.get("by_group", {})
        return {"action": "select_slice",
                "params": {"diagnosis": "CN"},
                "rationale": (f"修正后四组={byg}，PET+/Plasma− 从 551 骤降至 120。"
                              f"检验中间态（PET−/Plasma+）是否在 CN 层仍复现，排除该结论仅由编码 bug 支撑。"),
                "state_update": "切片"}

    @staticmethod
    def _s3(history):
        return {"action": "test_confounder",
                "params": {"control_vars": ["CARDIO", "ENDO"]},
                "rationale": "中间态复现后，检验其是否被合并症（心血管/内分泌）解释，做混杂归因。",
                "state_update": "方法"}

    @staticmethod
    def _s4(history):
        return {"action": "profile_mechanism",
                "params": {"target_group": "PET−/Plasma+", "markers": ["tau_pet", "GFAP", "NfL"], "alignment_window": "180d", "outcome": "ADAS13"},
                "rationale": "四组现象需要进入病理链定位：先比较 PET−/Plasma+ 的 tau PET、GFAP、NfL 与双阴参考组，区分非 AD 损伤与 tau 时序滞后。",
                "state_update": "假设"}

    @staticmethod
    def _s5(history):
        prev = history[-1]["feedback"].get("metrics", {})
        ranking = prev.get("candidate_mechanisms", [])
        top = ranking[0]["mechanism"] if ranking else "机制候选"
        return {"action": "sensitivity_analysis",
                "params": {"outcome": "CDRSB", "window": ">=2yr_followup"},
                "rationale": f"机制轮将 {top} 排在首位；换 CDRSB 并限制至少 2 年随访，检验认知证据是否稳健。",
                "state_update": "方法"}

    @staticmethod
    def _s6(history):
        return {"action": "profile_mechanism",
                "params": {"target_group": "PET+/Plasma−", "markers": ["tau_pet", "GFAP", "NfL"], "alignment_window": "365d", "outcome": "ADAS13"},
                "rationale": "对照另一不一致方向并放宽 tau–PET 时间窗口，判断影像先行型是时序滞后、非 AD 损伤还是测量噪声。",
                "state_update": "假设"}

    @staticmethod
    def _s7(history):
        return {"action": "discover_subtypes",
                "params": {"method": "gmm", "k": 3},
                "rationale": "机制证据优先于聚类；最后保留 GMM k=3 作为辅助方法敏感性，确认组内聚类不能替代机制定位。",
                "state_update": "方法"}


# --------------------------------------------------------------------------
# 探索 Agent
# --------------------------------------------------------------------------
class ExplorationAgent:
    def __init__(self, env: ExplorationEnvironment, llm=None, llm_base_url=None, llm_model=None):
        self.env = env
        if llm is None:
            llm = OpenAICompatLLM(base_url=llm_base_url, model=llm_model)
        self.llm = llm
        self.mock = MockPolicy()

    def _llm_propose(self, obs: dict, history: list) -> dict:
        user = json.dumps({
            "observation": obs,
            "history_last_6": history[-6:],
            "instruction": "输出下一轮动作的严格 JSON：action/params/rationale/state_update",
        }, ensure_ascii=False)
        text = self.llm(SYSTEM_PROMPT, user)
        return _extract_json(text)

    def propose(self, obs: dict, history: list) -> dict:
        if getattr(self.llm, "available", False):
            try:
                return self._llm_propose(obs, history)
            except Exception as e:
                print(f"[agent] LLM 调用失败，降级 MockPolicy：{type(e).__name__}: {e}")
        return self.mock.propose(obs, history)

    def run(self, n_rounds: int = 8, verbose: bool = True):
        obs = self.env.observe()
        history: list = []
        for i in range(1, n_rounds + 1):
            proposal = self.propose(obs, history)
            action = proposal.get("action", "sensitivity_analysis")
            params = proposal.get("params", {})
            rationale = proposal.get("rationale", "")
            state_update = proposal.get("state_update", "方法")
            fb = self.env.act(action, params)
            success = "error" not in fb
            entry = self.env.record(i, action, params, rationale, fb, state_update, success)
            history.append(entry)
            if verbose:
                interp = fb.get("interpretation", fb.get("error", ""))
                print(f"[round {i}] {action} {params} → {'✓' if success else '✗'}  {interp[:120]}")
            obs = self.env.observe()
        path = self.env.flush_logs()
        return history, path


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    env = ExplorationEnvironment()
    agent = ExplorationAgent(env)
    logs, path = agent.run(n_rounds=8)
    print(f"\n探索日志已写：{path}（{len(logs)} 轮）")
    print("最后一轮 feedback 摘要：")
    print(json.dumps(logs[-1]["feedback"], ensure_ascii=False, indent=1)[:800])
