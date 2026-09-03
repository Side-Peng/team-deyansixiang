# 影像–血浆生物标志物「不一致」亚型的 Agent 探索环境：跨切片复现性与认知轨迹关联

GOAI AI for Research 开放探索赛作品（初赛 + 复赛）。研究问题：**在血浆标志物与 Aβ PET 影像证据不一致的认知障碍人群中，是否存在跨人群稳定复现、且认知轨迹可区分的不一致亚型？**

## 复赛（半决赛）交付：三件套

复赛要求「最小可运行探索环境 + 探索日志 + 参照系设计」，对应文件：

| 交付物 | 文件 | 说明 |
|---|---|---|
| ① 最小可运行探索环境 | `scripts/environment.py` | `observe()` / `act()` / `record()` 统一接口；5 个动作（define_discordance / discover_subtypes / select_slice / test_confounder / sensitivity_analysis）；每轮统一返回 5 项反馈（Page 2 §2.2） |
| ② 探索日志 | `data/processed/exploration_log.jsonl` + `exploration_report.md` | 逐轮 round_log（round_id/action/params/rationale/state_update/feedback），失败轮保留 |
| ③ 参照系设计 | `scripts/reference.py` → `data/processed/reference_report.txt` | 平凡解 / 随机参照 / 无干预 / Pyun 四组 baseline 四层（Page 3 §3.2） |
| Agent（闭环） | `scripts/agent.py` | LLM 驱动（OpenAI 兼容，`LLM_API_KEY` 可配），无 key 自动降级确定性 MockPolicy |
| 一键复现 | `scripts/run_exploration.py` | 跑通三件套 + 生成人可读报告 |
| 可视化仪表盘 | `dashboard.html`（`scripts/make_dashboard.py` 生成） | 双击打开，离线可视化 8 轮闭环 + 参照系 + 修正后结论 |

运行（在 `goai-open-exploration/` 下，先装依赖，见下）：

```bash
python scripts/run_exploration.py                # 无 key：确定性 MockPolicy
LLM_API_KEY=sk-xxx python scripts/run_exploration.py   # LLM 驱动闭环
```

依赖（Windows 下按 CLAUDE.md 设 `PYTHONIOENCODING=utf-8`）：

```bash
python -m venv .venv && .venv/bin/pip install -r scripts/requirements.txt
# Windows 激活：.venv\Scripts\activate
```

### 复赛探索的诚实记录（重要）

在构建复赛探索环境、复现初赛口径时发现一处**口径 bug**：初赛 `preprocess.py` 用 `PTAU217 > 0.183` 生成 `PLASMA_STATUS`，pandas 中 `NaN > 阈值 = False`，导致**「未测血浆」的个体被错误编码为「血浆阴性」**。环境内建 `missing_policy` 参数（默认 `exclude` 正确口径、`missing_as_negative` 可复现初赛口径）；LLM Agent 的 8 轮自主探索从修正口径出发，独立确认了修正后结论。

- 修正口径（缺失排除）：不一致率 **18.4%**（248/1349，与初赛 `calibrate_report.txt` 标定一致）；四组 = 641 / 128 / 120 / 460。
- 初赛口径（缺失按阴性）：不一致率 31.6%；四组 = 1012 / 128 / 551 / 460。其中「PET+/Plasma− 551 人」里 **431 人（78%）实际是未测血浆**，真「血浆阴性」仅 120 人。

这一修正**不推翻核心发现**：中间态（PET−/Plasma+，128 人）不受影响（该组要求血浆阳性，天然排除缺失），其跨切片复现、相对双阴的轨迹差异（修正后 p=0.0325）与混杂校正后保留均成立。但它意味着初赛对「PET+/Plasma− 影像先行型」的轨迹/tau 解读需在排除缺失后复核。详见探索日志 Round 1–2 与 `four_page_guide.md`。

**tau 复核结果（`scripts/tau_analysis_corrected.py` → `data/processed/tau_report_corrected.txt`）**：修正口径下，PET+/Plasma− 组 tau 可评估覆盖从初赛的 26.7%（147/551）升到 **80.8%**（97/120，样本更干净），但 tau 负担明显下调——meta-temporal SUVR 从 **1.230 → 1.195**，≥1.28 阳性率从 36.7% → 18.6%。修正后两个不一致方向的 tau 负担**几乎对称**（PET+/Plasma− 1.195 ≈ PET−/Plasma+ 1.191），都仅略高于双阴（1.169，p≈0.002–0.005）且远低于双阳（1.387）。**结论修正**：初赛「PET+/Plasma− = tau 已启动的真实病理早期」的解读在排除缺失后**不再成立**，两个不一致方向在 tau 维度上趋同（都接近「无显著 tau 启动」），初赛的「不对称/非镜像」叙事需下调。

## 机制版管线（2026-08-25）

完整研究版一键运行：

```bash
python scripts/run_exploration.py 8
```

新增产物：

- `data/processed/raw_validation.json`：原始表 schema 检查
- `data/processed/mechanism_report.json`：Aβ PET → plasma p-tau217 → tau PET → GFAP/NfL → cognition 机制证据
- `data/processed/validation_report.json`：按 PET site 的 patient-level held-out validation
- `dashboard.html`：机制排序、时间对齐、缺失率、验证结果和 Agent 探索时间线

机制动作 `profile_mechanism` 不训练 protein foundation model，而是使用真实的 GFAP/NfL 蛋白标志物测量，比较 `temporal_lag`、`non_ad_injury` 和 `measurement_noise` 三类候选机制。所有结果均为假设生成，报告同时保留效应量、95% CI、缺失率、时间间隔和 Benjamini-Hochberg q 值。

API key 只通过 `LLM_API_KEY`、`LLM_BASE_URL`、`LLM_MODEL` 配置；`llm_config.json` 不应保存真实 secret。ADNI raw data 和个体级宽表不入库。

## 本地探索台（Web 界面）

```bash
python webapp/server.py        # → http://127.0.0.1:8765
```

零第三方依赖（标准库后端 + 原生前端）。四个面板：**总览**（四组结构/轨迹/tau 负担实时统计 + 探索时间线）、**LLM 探索**（OpenAI 兼容接口驱动探索闭环，无 key 自动降级 MockPolicy）、**数据**（拖拽上传 CSV 替换管线表 + 一键重跑管线）、**报告**（浏览全部聚合报告）。数据源默认 ADNIMERGE2 总和数据（`data/raw/merged/`，由 `scripts/export_adnimerge2.R` 导出），上传替换后需在「数据」页重跑管线生效。

## 项目结构

```
├── four_page_guide.md                # 工作区母稿（证据台账/探索记录/未决问题）
├── four_page_guide_submission.md     # 提交版四页稿（源文件）
├── four_page_guide_submission.docx   # 提交版 Word（四页自动分页）
├── scripts/                          # 完整可复现管线（7 个脚本 + 配置）
│   ├── config.py                     # 阈值常量 + 字段映射（集中配置）
│   ├── preprocess.py                 # 清洗 + 基线定义 + 宽表构建
│   ├── calibrate.py                  # Youden 阈值标定（含逆向标志物方向处理）
│   ├── define_discordance.py         # 四组切分 + E1 定性核对
│   ├── trajectory.py                 # 纵向认知轨迹（KW/校正回归/置换检验）
│   ├── slice_analysis.py             # 合并症混杂 + 诊断分层复现
│   ├── cluster_pm.py                 # P−/P+ 组内聚类（双亚群验证）
│   ├── tau_analysis.py               # tau PET 负担机制区分
│   └── md2docx.py                    # 文档转 Word 工具（无依赖）
└── data/processed/                   # 全部结果报告（聚合统计）
    ├── calibrate_report.txt          # 阈值标定
    ├── group_stats.csv               # 四组描述统计
    ├── trajectory_report.txt         # 轨迹
    ├── slice_report.txt              # 合并症+分层
    ├── cluster_report_v2.txt         # 聚类
    └── tau_report.txt                # tau 负担
```

## 运行（三步）

1. 按 `scripts/config.py` 中 FILES 映射，将 ADNI 下载的表置于 `data/raw/`（临床表 + UPENN 血浆 + UC Berkeley amyloid/tau PET 处理表）
2. 安装依赖：`pip install pandas numpy scipy scikit-learn`
3. 依次运行（每步输出落盘，可逐轮复现）：

```bash
python scripts/preprocess.py            # → subjects_wide.csv
python scripts/calibrate.py            # → 阈值（已标定，复核后回填 config）
python scripts/define_discordance.py   # → group_stats.csv（E1 核对 3/3）
python scripts/trajectory.py           # → trajectory_report.txt
python scripts/slice_analysis.py       # → slice_report.txt
python scripts/cluster_pm.py           # → cluster_report_v2.txt
python scripts/tau_analysis.py         # → tau_report.txt
```

## 五轮试跑结果（2026-08-02/03，ADNI 全量已下载）

1. **标定+四组**：p-tau217 vs Aβ PET AUC 0.885，Youden 阈值 0.183 pg/mL，不一致率 18.4%（同 E4 量级），四组结构复现 E1
2. **纵向轨迹**：PET+/Plasma−（551 人）ΔADAS13/年 +1.31 接近双阳（+1.66）；PET−/Plasma+（128 人）温和中间态 +0.46（vs 双阴 p=0.006，vs 双阳 p=1.5e-11）；置换 p=0.000
3. **合并症+分层**：PET+/Plasma− 合并症最高（心血管 49%/内分泌 36%，p=3e-20）；CN 层中间态复现
4. **聚类**：合并症亚群可切出；两亚群轨迹分化不显著（p=0.68，部分负结果，如实报告）
5. **tau 负担**：PET+/Plasma− SUVR 1.230 显著高于双阴（p=4.5e-12）低于双阳（p=7.1e-12）→ 真实病理早期；PET−/Plasma+ SUVR ≈ 双阴 → 无 tau 启动（非 AD 途径为主）

核心发现：**不一致的两个方向不是镜像**——影像先行型（PET+/Plasma−）= 真实病理早期、血浆滞后；血浆孤立型（PET−/Plasma+）= 以非 AD 途径为主。

## 数据与合规

- ADNI 数据受 Data Use Agreement 约束，**原始数据与个体级宽表（subjects_wide.csv）不随本包分发**；本包仅含聚合统计报告与复现脚本
- 本地 150 例数据按医院伦理与脱敏要求处理，仅用于外部验证（下一阶段）
- 过程中发现并修复一个标签命名 bug（GROUP 1/2 组名对调），修复过程见 four_page_guide.md 诚实记录，全部结论为修正后重跑结果
