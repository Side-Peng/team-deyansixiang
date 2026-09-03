# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目定位

GOAI AI for Research 开放探索赛初赛作品。研究问题：血浆 p-tau217 与 Aβ PET 阳性状态不一致的认知障碍人群中，是否存在跨切片稳定复现、与认知轨迹/tau 病理相关的不一致亚型。数据主体为 ADNI（已获批，受 DUA 约束），本地 150 例仅外部验证。

母稿 `four_page_guide.md`（证据台账/探索记录/未决问题）与提交版 `four_page_guide_submission.md/.docx` 已**移出仓库**，在 `E:\2026GOAI\` 下；仓库内保留的探索记录是 `semifinal_exploration_record.md`（复赛口径 bug 发现与修正的证据链）。

## 环境注意

- **Windows + Git Bash**：脚本必须设 `export PYTHONIOENCODING=utf-8`，否则 GBK 编码会在含中文的 print/写文件时崩溃。
- 运行脚本时工作目录应在 `goai-open-exploration/`，`import config` 依赖 `scripts/` 在 path 上（各脚本从 `scripts/` 同级运行：`python scripts/preprocess.py`）。
- ADNI 原始 CSV 位于 `data/raw/`，**不入库**（`.gitignore` 已排除，DUA 约束）；`data/processed/subjects_wide.csv` 含 RID 个体级数据，同样不入库，分发仅限聚合报告（`*_report.txt`、`group_stats.csv`）。

## 数据源：ADNIMERGE2 总和数据（2026-09-02 切换，替代 IDA 单表 CSV）

- 原数据：`data/raw/*.csv`（IDA 逐表下载，2026-08-02 版）。现数据：`ADNI_MERGE/ADNIMERGE2.tar.gz`（R 包，2026-01-05 版，含 ADNI1/GO/2/3/4）。
- 转换：`scripts/export_adnimerge2.R` 用 R 4.4.3（`S:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe`，base R 即可）把管线需要的 10 张表从 `ADNIMERGE2/data/*.rda` 导出为 `data/raw/merged/<TABLE>.csv`；`config.FILES` 已指向 `merged/`。列名与老 IDA 表逐一核对一致，`FIELDS` 未动。
- **切换时踩过的坑（均已修）**：① ADNIMERGE2 把值重编码为字符串（PTGENDER 'Male'/'Female'、DIAGNOSIS 'CN'/'MCI'/'Dementia'、MH4CARD 'Yes'/'No'）→ config 三个 normalize 函数统一；② `PTDOBYY` 是 float 年份，`pd.to_datetime(1931.0)` 会被当作纳秒时间戳 → 1970 年，AGE 全错，改为年份数值直接相减；③ `calibrate.py` 两个标志物共用 dropna 样本，AB42/40 缺失模式与 pT217 不完全重合 → 每标志物单独 dropna；④ `cluster_pm.make_comorbidity` 的 CKD 未二值化（INITHEALTH 命中行数直接当 0/1，阳性率可 >100%）→ 已修 `>0`。
- **新数据四组基线**（修正口径 exclude，校准样本 1134）：541 / 103 / 99 / 391，不一致率 17.8%（Youden 0.183，AUC 0.889，与老数据标定一致）。老数据是 641 / 128 / 120 / 460（1349）。关键结论方向不变：MCI 层严格中间态复现（p=0.021）、PET+/Plasma− 修正口径 tau SUVR 1.195 ≈ 双阴水平、轨迹 +0.35 接近双阴。
- `ADNI_MERGE/` 与 `data/raw/` 均已 gitignore（DUA）。

## 管线结构

8 个脚本构成线性管线，全部配置集中在 `scripts/config.py`（`FILES` 文件名映射、`FIELDS` 列名映射、`THRESHOLDS` 阈值）。每步读上一步产物或原始表，输出落盘到 `data/processed/`，无统一 CLI 参数——改实验改 `config.py`。

```
preprocess.py        原始 CSV → subjects_wide.csv（基线宽表）
calibrate.py         Youden 阈值标定（AUC/四格表）→ calibrate_report.txt/json
define_discordance.py  四组切分 + E1 定性核对 → group_stats.csv
trajectory.py        纵向年化变化率 + KW/MW/校正回归/置换检验
slice_analysis.py    合并症混杂 + CN/MCI/AD 诊断分层复现
cluster_pm.py        PET−/Plasma+ 组内 KMeans 双亚群（含 INITHEALTH 关键词合并症）
tau_analysis.py      tau PET meta-temporal SUVR 机制区分（FTP 示踪剂分层）
md2docx.py           提交版 Markdown → .docx（纯 OOXML zip，无第三方依赖）
```

无 lint 框架；测试在 `tests/`（pytest，`pyproject.toml` 配置了 `pythonpath=scripts`）。验证方式：关键分析逐脚本运行并人工核对 `data/processed/*_report.txt` 中的数字。

## 关键约定与易踩坑点（跨文件才看得出来）

- **基线定义**：每张表按 RID 取**日期最早行**（`earliest()` in preprocess.py），**不依赖 VISCODE='bl'**——ADNI 各表 VISCODE 编码不统一（CDR/MMSE 用 'sc'，PTDEMOG 每人多行）。新增表必须走同一逻辑。
- **GROUP 编码**：`GROUP = PET_STATUS*2 + PLASMA_STATUS`，即 **GROUP=1 是 PET−/Plasma+，GROUP=2 是 PET+/Plasma−**。曾因三个脚本的命名表把 1/2 对调导致结论错误（2026-08-03 已修复，母稿有诚实记录）。新增组间比较时务必按此数值映射，不要凭直觉命名。
- **逆向标志物**：Aβ42/40 比值**越低越病理**，`calibrate.py` 的 `youden()` 用 `direction=-1`（四格表用 `s < thr` 而非 `>`）。roc_auc_score 默认"分数越高越正类"，逆向标志物不处理会得到 AUC<0.5。
- **ADNI 缺失编码**：`-4.0`、`≤0` 在 CLEAN_RULES 中清洗为 NA。**编码规范统一在 config**：`PTGENDER`/`DIAGNOSIS`/MEDHIST 标记在 ADNIMERGE2 是字符串（'Male'/'Female'、'CN'/'MCI'/'Dementia'、'Yes'/'No'），老 IDA 是数字（1/2、1/2/3、0/1）；`config.normalize_gender/normalize_dx/normalize_yes_no` 统一为数字编码，下游回归按数字消费（GENDER map `{1:1, 2:0}`）。新数据源表列名与老表一致，但值编码不同，改数据源时先查这四处。
- **tau 分析**：只取 `TRACER=='FTP'` 且 `qc_flag∈{1,2}`，每 RID 最早扫描——表自带警告 "DO NOT COMPARE SUVRs ACROSS TRACERS"，MK6240/PI2620 不能与 FTP 混算。tau 阳性阈值（1.22/1.28/1.36）是文献常用 cutoff，**非 ADNI 官方判定**，报告中必须注明。
- **合并症双来源**：MEDHIST 系统标记（MH4CARD/MH9ENDO，仅 ADNI1/GO/2）∪ INITHEALTH 文本关键词（含否定短语排除，见 `cluster_pm.py` 的 `make_comorbidity()`）。
- **置换检验**：`trajectory.py` 用 1000 次打乱组标签重算 KW H；种子在脚本内 `np.random.default_rng(2026)`，与 config 的 SEED=42 不同（历史原因，改时注意）。
- **`compute_plasma_status` 有两份**：`environment.py` 与 `data_layer.py` 各实现一份（`corrected_reports.py` 用 environment 版，`mechanism_analysis.py`/`validation.py` 用 data_layer 版），口径逻辑改动必须两处同步。
- **GFAP_Q / NfL_Q 是 `data_layer.build_aligned_cohort` 硬编码的 plasma 表列名**（未进 `config.FILES`），换数据源/换表时此处先查是否仍存在（缺失 → 整列 NA，机制分析悄悄变空）。
- **机制层要 `PET_SITEID`**：`validation.py` 依赖 pet 表的 `SITEID` 列（preprocess 主管线不导出该列），数据源切换后先确认。
- **md2docx.py**：自写轻量 Markdown 子集 → OOXML（python-docx/pandoc 在本机不可用，Word COM 在非交互会话报 0x80080005）。h2 自动 `pageBreakBefore`（四页结构），表格首行加粗。提交 docx 改后用 `python scripts/md2docx.py four_page_guide_submission.md four_page_guide_submission.docx` 重新生成。

## 依赖

`pip install pandas numpy scipy scikit-learn`（依赖清单见 `pyproject.toml`，锁定版 `requirements-lock.txt`，旧 `scripts/requirements.txt` 仍在）。FreeSurfer 仅在做 MRI 影像预处理时需要，当前管线未启用（MRI 体积作协变量为设计预留）。

## 复赛探索环境（2026-08-25 新增，区别于初赛静态管线）

复赛新增 5 个脚本，构成「Agent 探索环境」层，**不改动**初赛 8 个脚本（仅 `cluster_pm.py` 做了一处 pandas 3 兼容修复，见下）：

- `environment.py` —— `ExplorationEnvironment`：`observe()` / `act(action, params)` / `record(...)`；5 个动作对应 Page 2 §2.2 的行动，每个动作返回统一 5 项反馈（①不一致组规模 ②跨切片复现系数 ③轨迹分离效应量+CI+p ④混杂校正后变化 ⑤置换参照）。复用初赛 `trajectory.per_patient_change`、`cluster_pm.make_comorbidity`、`calibrate.youden`。
- `reference.py` —— `ReferenceFrame`：Page 3 §3.2 四层参照（平凡解/随机/无干预/Pyun baseline）。
- `agent.py` —— `ExplorationAgent`：LLM 驱动闭环（`OpenAICompatLLM`，标准库 urllib，读 `LLM_API_KEY`/`LLM_BASE_URL`/`LLM_MODEL`）；无 key 时 `MockPolicy`（确定性、读真实反馈、可复现）。
- `run_exploration.py` —— 一键产出：`raw_validation.json`（原始表 schema 校验，失败即退出）→ `exploration_log.jsonl` + `exploration_report.md` + `reference_report.txt` → `mechanism_report.json` + `validation_report.json`。
- `make_dashboard.py` —— 生成自包含静态仪表盘 `dashboard.html`（双击打开、离线可用，数据内嵌自 `exploration_log.jsonl` + `reference_report.json` + 机制报告）。跑完 `run_exploration.py` 后重跑本脚本即可刷新。

### 机制画像层（2026-08-25 增量）

- `data_layer.py` —— canonical data layer：`build_baseline()`（同行基线、**保留缺失**，供修正口径）、`compute_plasma_status`/`add_canonical_group`（`missing_policy="exclude"` 为修正口径）、`build_aligned_cohort(window∈{90,180,365}d)`（血浆/GFAP/NfL 与 tau PET 用 merge_asof nearest 对齐到 PET 扫描日，`*_TIME_GAP_DAYS` 记录缺口）、`summarize_missingness`、`validate_raw_inputs`。
- `mechanism_analysis.py` —— 机制画像：temporal lag / non-AD injury（GFAP+NfL）/ measurement noise 三类候选机制证据 + 评分排序，HuberRegressor 校正轨迹差异，Benjamini-Hochberg 校正 → `mechanism_report.json`。
- `validation.py` —— internal validation：按 `PET_SITEID` 每 5 站点留一，RID 不跨 split，固定阈值，比较 train/valid 的机制排序 top-2 重合率 → `validation_report.json`。
- `corrected_reports.py` —— 修正口径（缺失血浆排除）重跑 trajectory / slice / cluster → `*_report_corrected.txt` 三份。
- `pyproject.toml` + `tests/`：pytest 配置（`pythonpath=scripts`），`test_data_layer.py` / `test_environment.py`，跑 `pytest` 验证；`requirements-lock.txt` 为锁定版依赖。
- `webapp/start.bat`：Windows 双击启动器（先探测 8765 端口防重复启动，再开浏览器）。

关键结果：机制排序两方向一致（non_ad_injury > temporal_lag > measurement_noise），site-held-out top-2 重合率两方向均 1.0；但 180d 对齐 cohort 仅 717 人（P−/P+ 41 / P+/P− 62），排序只是假设生成，报告含 CI/缺失率/q 值。

运行：本仓库在 WSL 下用 `.venv`（`python3 -m venv .venv && .venv/bin/pip install -r scripts/requirements.txt`），入口 `.venv/bin/python scripts/run_exploration.py`。Windows 原生环境同样适用（激活 `.venv\Scripts\activate`，记得 `PYTHONIOENCODING=utf-8`）。

## 本地探索台 Web 界面（2026-09-02 新增，零第三方依赖）

`webapp/`（`server.py` + `index.html` + `app.js` + `style.css`）——标准库 ThreadingHTTPServer 后端 + 原生 JS/Canvas 前端，`python webapp/server.py` → http://127.0.0.1:8765。四个面板：

- **总览**：四组结构 / 轨迹 / tau 负担（pandas 现算 + 60s 缓存，任务完成自动失效）、探索时间线
- **LLM 探索**：OpenAI 兼容接口驱动 `agent.ExplorationAgent`（无 key 降级 MockPolicy）；「记住配置」写 `llm_config.json`
- **数据**：上传 CSV 覆盖 `data/raw/merged/` 同名表（白名单 = FILES 的 basename），⚡快速重跑（preprocess→四组）/ ▶完整重跑管线
- **报告**：浏览 `data/processed/*.txt`

易踩坑（已修）：① `start_task` 持锁时调 `_set_task`（内部也拿锁）→ 非可重入 Lock 死锁，必须 `RLock`；② `compute_stats` 里 `group` 是 np.int64，json.dumps 前要 `int()`；③ 上传白名单要用 `Path(f).name`（FILES 值带 `merged/` 前缀）；④ 多实例绑同端口（Windows SO_REUSEADDR），重启前 `netstat -ano | grep 8765` 杀干净。

### ⚠️ 关键口径 bug（复赛探索 Round 1–2 发现，影响初赛结论解读）

初赛 `preprocess.py` 生成 `PLASMA_STATUS = (PTAU217_bl > 0.183)`，pandas 里 `NaN > 阈值 = False`，**未测血浆的个体被编码为「血浆阴性」**：

- 正确口径（`environment.compute_plasma_status(..., missing_policy="exclude")`）：不一致率 18.4%，四组 641/128/120/460（与 `calibrate_report.txt` 一致）。
- 初赛口径（`missing_policy="missing_as_negative"`）：不一致率 31.6%，四组 1012/128/551/460；「PET+/Plasma− 551」中 431 人（78%）是未测血浆。

结论：中间态（PET−/Plasma+ 128 人）与双阳（460）不受影响（需血浆阳性）；但初赛「PET+/Plasma− 影像先行型」的轨迹（+1.31）与 tau 负担（SUVR 1.230）是在混入 78% 未测血浆的组上算的。**修正口径复核已落盘**：`corrected_reports.py` 产出 `trajectory_report_corrected.txt` / `slice_report_corrected.txt` / `cluster_report_corrected.txt`，`tau_analysis_corrected.py` 产出 `tau_report_corrected.txt`。修正后 PET+/Plasma− 轨迹中位 +0.35 ≈ 双阴 +0.23（远离双阳 +1.66），MCI 层严格中间态仍复现（p=0.021）。写复赛报告时务必引用修正口径。

**tau 复核（2026-08-25，`tau_analysis_corrected.py`）**：修正口径下 PET+/Plasma− 组 tau 可评估 97/120（80.8%，初赛口径 147/551=26.7%），SUVR 1.195（初赛 1.230）、≥1.28 阳性率 18.6%（初赛 36.7%）。两个不一致方向 tau 负担趋同（1.195 ≈ 1.191，都≈双阴 1.169、远低于双阳 1.387）→ **初赛「PET+/Plasma− = tau 已启动的早期病理」与「两方向不对称」结论被修正为「两方向在 tau 维度趋同、均无显著 tau 启动」**。产出 `data/processed/tau_report_corrected.txt`。

### pandas 3 兼容修复

`cluster_pm.make_comorbidity` 的 INITHEALTH 拼接改为先 `fillna("")` 再 `astype(str)`：pandas ≥3 下 `astype(str)` 对 NaN 保留 NaN（而非旧版 "nan"），会导致 `k in s` 对 float 报错；`fillna("")` 与旧版 "nan" 均不命中关键词，结果等价。
