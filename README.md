# 影像–血浆生物标志物「不一致」亚型的 Agent 探索环境

> 面向认知障碍人群亚型研究的可复现探索基础设施。聚焦血浆 p-tau217 与 Aβ PET 阳性状态不一致的人群,以「统一探索接口 + LLM 自主 Agent + 多层级参照系」为骨架,对不一致亚型的稳定性、认知轨迹与 tau 病理特征进行端到端、可审计的系统性探索。

## 研究问题

在血浆标志物与 Aβ PET 影像证据**不一致**的认知障碍人群中,是否存在跨人群稳定复现、且认知轨迹可区分的不一致亚型?

## 核心特性

- **统一探索接口** —— `observe()` / `act()` / `record()` 三层抽象,5 个探索动作(`define_discordance` / `discover_subtypes` / `select_slice` / `test_confounder` / `sensitivity_analysis`),每轮动作统一返回五项结构化反馈,所有探索路径均可回放与审计
- **LLM 驱动的自主 Agent** —— OpenAI 兼容接口,自动编排探索动作、记录推理与状态更新;未配置密钥时自动降级为确定性策略,保证零依赖可复现
- **四层参照系设计** —— 平凡解 / 随机参照 / 无干预 / 文献基线,让每一组结论都有锚点,而非孤立数字
- **机制画像与站点级验证** —— 三类候选机制(时间滞后 / 非 AD 损伤 / 测量噪声)证据排序,按 PET 站点做 held-out 验证,报告含效应量、置信区间与 FDR 校正
- **自包含可视化仪表盘** —— 单文件 `dashboard.html`,离线双击打开,汇总探索时间线、四组结构、轨迹与 tau 负担
- **隐私合规** —— 数据层与报告层分离,仅分发聚合统计与各阶段报告,不含任何个体级数据

## 快速开始

```bash
# 1. 安装依赖 (Python ≥ 3.10)
pip install -r scripts/requirements.txt
# 或: pip install -e .

# 2. 一键运行探索(确定性模式)
python scripts/run_exploration.py

# 3. LLM 驱动模式(可选)
export LLM_API_KEY=sk-xxx
export LLM_BASE_URL=https://api.openai.com/v1
export LLM_MODEL=gpt-4o
python scripts/run_exploration.py 8
```

Windows 使用 Git Bash 时请先 `export PYTHONIOENCODING=utf-8`。

## 使用说明

### 命令行

`python scripts/run_exploration.py [探索轮数]`

一键产出:

| 产物 | 说明 |
|---|---|
| `exploration_log.jsonl` | 逐轮探索日志(动作 / 参数 / 推理 / 状态更新 / 反馈),失败轮保留 |
| `exploration_report.md` | 人可读探索月报 |
| `reference_report.txt` | 参照系对比 |
| `mechanism_report.json` | 机制画像证据排序 |
| `validation_report.json` | 站点级 held-out 验证 |
| `raw_validation.json` | 原始表 schema 校验结果 |

### Web 探索台

```bash
python webapp/server.py    # → http://127.0.0.1:8765
```

零第三方依赖(标准库后端 + 原生前端),四个面板:

- **总览** —— 四组结构 / 轨迹 / tau 负担实时统计 + 探索时间线
- **LLM 探索** —— OpenAI 兼容接口驱动闭环(无密钥自动降级)
- **数据** —— 拖拽上传替换管线数据表,一键重跑
- **报告** —— 浏览全部聚合报告

### 运行测试

```bash
pytest
```

## 配置

| 环境变量 | 说明 |
|---|---|
| `LLM_API_KEY` | LLM 密钥,通过环境变量注入,不落盘 |
| `LLM_BASE_URL` | OpenAI 兼容端点 |
| `LLM_MODEL` | 模型名 |

参考模板见 `llm_config.example.json`。真实密钥**不要**写入仓库内的任何配置文件。

## 项目结构

```
├── scripts/
│   ├── environment.py          # 探索环境:observe()/act()/record() + 5 个动作
│   ├── agent.py                # LLM 驱动闭环 Agent(无 key 降级 MockPolicy)
│   ├── reference.py            # 四层参照系
│   ├── run_exploration.py      # 一键入口
│   ├── make_dashboard.py       # 生成 dashboard.html
│   ├── data_layer.py           # 规范化数据层(基线/口径/对齐)
│   ├── mechanism_analysis.py   # 机制画像
│   ├── validation.py           # 站点级验证
│   ├── corrected_reports.py    # 修正口径复核报告
│   ├── preprocess.py …         # 原始数据 → 四组切分的核心管线
│   └── config.py               # 阈值与字段映射集中配置
├── webapp/                     # 本地探索台(标准库 Web)
├── tests/                      # pytest 测试
├── dashboard.html              # 自包含仪表盘(生成产物)
└── data/processed/             # 全部聚合报告与日志
```

## 数据与合规

数据源为 ADNI(阿尔茨海默病神经影像计划),受**数据使用协议(DUA)**约束:

- 原始数据与个体级宽表**不随本仓库分发**
- 本仓库仅含聚合统计、报告文本与复现脚本
- 替换自有数据源时,将按 `scripts/config.py` 中 `FILES` 映射命名的 CSV 置于 `data/raw/`,即可复用完整管线
