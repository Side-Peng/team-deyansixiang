# 探索报告（GOAI AI4S 探索赛道 · 复赛）

> 由 `scripts/run_exploration.py` 自动生成；逐轮可复现，日志见 `exploration_log.jsonl`。

| 轮 | 动作 | 参数 | 状态更新 | 关键结果 |
|---|---|---|---|---|
| 1 | define_discordance | `{"assay_pair": "ptau217_vs_pet", "missing_policy": "missing_as_negative"}` | 定义 | 不一致率 32.8%；四组 {'PET−/Plasma−': 1018, 'PET−/Plasma+': 103, 'PET+/Plasma−': 585, 'PET+/Plasma+': 391} |
| 2 | define_discordance | `{"assay_pair": "ptau217_vs_pet", "missing_policy": "exclude"}` | 定义 | 不一致率 17.8%；四组 {'PET−/Plasma−': 541, 'PET−/Plasma+': 103, 'PET+/Plasma−': 99, 'PET+/Plasma+': 391} |
| 3 | select_slice | `{"diagnosis": "CN"}` | 切片 | 中间态=成立（p_int_vs_neg=0.232） |
| 4 | test_confounder | `{"control_vars": ["CARDIO", "ENDO"]}` | 方法 | Δ系数=-0.03404 |
| 5 | profile_mechanism | `{"target_group": "PET−/Plasma+", "markers": ["tau_pet", "GFAP", "NfL"], "alignment_window": "180d", "outcome": "ADAS13"}` | 假设 | 机制排序=non_ad_injury, temporal_lag, measurement_noise |
| 6 | sensitivity_analysis | `{"outcome": "CDRSB", "window": ">=2yr_followup"}` | 方法 | 复现系数=1 |
| 7 | profile_mechanism | `{"target_group": "PET+/Plasma−", "markers": ["tau_pet", "GFAP", "NfL"], "alignment_window": "365d", "outcome": "ADAS13"}` | 假设 | 机制排序=non_ad_injury, temporal_lag, measurement_noise |
| 8 | discover_subtypes | `{"method": "gmm", "k": 3}` | 方法 | 轮廓=0.41；轨迹 p=0.298 |

## 每轮 rationale

- **Round 1**（define_discordance，定义）：先按初赛管线口径（缺失血浆按阴性）复现四组基线，确认与已提交结果一致。
- **Round 2**（define_discordance，定义）：初赛口径不一致率=32.8%，显著高于标定口径 18.4% 与 E4 量级；marker_missing_by_pet={0: 477, 1: 486} 提示大量「未测血浆」被编码为阴性。改 missing_policy=exclude 复核定义，验证是否为编码 bug。
- **Round 3**（select_slice，切片）：修正后四组={'PET−/Plasma−': 541, 'PET−/Plasma+': 103, 'PET+/Plasma−': 99, 'PET+/Plasma+': 391}，PET+/Plasma− 从 551 骤降至 120。检验中间态（PET−/Plasma+）是否在 CN 层仍复现，排除该结论仅由编码 bug 支撑。
- **Round 4**（test_confounder，方法）：中间态复现后，检验其是否被合并症（心血管/内分泌）解释，做混杂归因。
- **Round 5**（profile_mechanism，假设）：四组现象需要进入病理链定位：先比较 PET−/Plasma+ 的 tau PET、GFAP、NfL 与双阴参考组，区分非 AD 损伤与 tau 时序滞后。
- **Round 6**（sensitivity_analysis，方法）：机制轮将 non_ad_injury 排在首位；换 CDRSB 并限制至少 2 年随访，检验认知证据是否稳健。
- **Round 7**（profile_mechanism，假设）：对照另一不一致方向并放宽 tau–PET 时间窗口，判断影像先行型是时序滞后、非 AD 损伤还是测量噪声。
- **Round 8**（discover_subtypes，方法）：机制证据优先于聚类；最后保留 GMM k=3 作为辅助方法敏感性，确认组内聚类不能替代机制定位。
