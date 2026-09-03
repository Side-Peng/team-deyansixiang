# -*- coding: utf-8 -*-
"""make_dashboard.py — 生成自包含静态 HTML 仪表盘（双击即可打开，离线可用）

读取 data/processed/exploration_log.jsonl + reference_report.json，把数据内嵌进
一个单文件 HTML（无需服务器、无需第三方库），输出到仓库根目录 dashboard.html。

用法：python scripts/make_dashboard.py
再跑 run_exploration.py / corrected_reports.py 后重跑本脚本即可刷新。
"""
import json
from pathlib import Path

import config as C

ROOT = Path(__file__).resolve().parents[1]


def load_log():
    logs = []
    with open(C.PROC_DIR / "exploration_log.jsonl", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                logs.append(json.loads(line))
    return logs


def load_ref():
    path = C.PROC_DIR / "reference_report.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def load_json(name):
    path = C.PROC_DIR / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


HTML = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AI4S 探索环境仪表盘 · GOAI 复赛</title>
<style>
:root{
  --ink:#06152e; --muted:#5b6b84; --line:#d8e6fb; --blue:#0b5cff; --green:#22b85a;
  --orange:#ff6b2a; --bg:#f7fcff; --card:#ffffff;
}
*{box-sizing:border-box}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",sans-serif;
  background:radial-gradient(circle at 82% 0,rgba(33,214,255,.16),transparent 30%),
  radial-gradient(circle at 10% 60%,rgba(68,214,154,.10),transparent 30%),var(--bg);
  color:var(--ink);line-height:1.6}
.wrap{max-width:1060px;margin:0 auto;padding:28px 20px 80px}
header.hero{background:linear-gradient(112deg,#12318d,#1743a2 30%,#0f5d92 57%,#061744);color:#fff;
  border-radius:20px;padding:30px 32px;margin-bottom:22px;position:relative;overflow:hidden}
header.hero h1{margin:0 0 8px;font-size:24px;font-weight:800;letter-spacing:-.01em}
header.hero p{margin:2px 0;opacity:.92;font-size:14px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-top:18px}
.kpi{background:rgba(255,255,255,.12);border:1px solid rgba(255,255,255,.22);border-radius:14px;padding:12px 14px}
.kpi b{display:block;font-size:22px;font-weight:800}
.kpi span{font-size:12px;opacity:.85}
h2{font-size:18px;margin:34px 0 12px;font-weight:800}
.card{background:var(--card);border:1px solid var(--line);border-radius:16px;padding:18px 20px;
  box-shadow:0 10px 30px rgba(13,118,255,.06)}
.timeline{position:relative;padding-left:26px}
.timeline:before{content:"";position:absolute;left:9px;top:6px;bottom:6px;width:2px;background:var(--line)}
.node{position:absolute;left:-26px;top:20px;width:18px;height:18px;border-radius:50%;background:#fff;
  border:3px solid var(--blue)}
.node.fix{border-color:var(--orange);background:var(--orange)}
.node.neg{border-color:#9aa7bd;background:#9aa7bd}
.round{position:relative;margin-bottom:14px}
.round .head{display:flex;align-items:center;flex-wrap:wrap;gap:8px}
.badge{font-size:12px;font-weight:700;padding:2px 10px;border-radius:999px;background:#eef3ff;color:var(--blue)}
.badge.fix{background:#fff0e8;color:var(--orange)}
.badge.neg{background:#eef1f6;color:#5b6b84}
.tag{font-size:12px;color:var(--muted);border:1px solid var(--line);border-radius:8px;padding:1px 8px}
.rationale{font-size:13.5px;color:var(--muted);margin:8px 0 4px;font-style:italic}
.key{font-size:13.5px;font-weight:700}
details{margin-top:8px;font-size:12px;color:var(--muted)}
details pre{background:#0f2346;color:#dbeafe;border-radius:10px;padding:12px;overflow:auto;
  max-height:320px;font-size:11.5px;line-height:1.5}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:12px}
.ref-card{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:14px 16px}
.ref-card h3{margin:0 0 6px;font-size:14px}
.ref-card .v{font-size:20px;font-weight:800;color:var(--blue)}
.ref-card p{margin:4px 0;font-size:12.5px;color:var(--muted)}
.callout{border-left:4px solid var(--orange);background:#fff8f3;border-radius:0 12px 12px 0;
  padding:14px 18px;margin-top:14px;font-size:13.5px}
.callout b{color:var(--orange)}
.foot{color:var(--muted);font-size:12px;margin-top:30px;text-align:center}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:#eef3ff;padding:1px 5px;border-radius:5px}
</style>
</head>
<body>
<div class="wrap">
  <header class="hero">
    <h1>影像–血浆「不一致」亚型 · Agent 探索环境仪表盘</h1>
    <p>GOAI AI for Research 探索赛道 · 复赛 | 修正口径（缺失血浆排除）</p>
    <div class="kpis" id="kpis"></div>
  </header>

  <h2>探索时间线（8 轮闭环）</h2>
  <div class="card"><div class="timeline" id="timeline"></div></div>

  <h2>参照系（四层，来自 reference_report.json）</h2>
  <div class="grid" id="ref"></div>

  <h2>修正后核心结论（取代初赛三处）</h2>
  <div class="card" id="findings"></div>

  <h2>病理链机制定位</h2>
  <div class="grid" id="mechanisms"></div>
  <div class="card" id="validation"></div>

  <div class="foot">由 scripts/make_dashboard.py 生成 · 数据源 exploration_log.jsonl / reference_report.json / mechanism_report.json / validation_report.json · 逐轮可复现</div>
</div>

<script>
const DATA = __DATA__;
const logs = DATA.logs || [];
const ref  = DATA.ref  || {};
const mechanism = DATA.mechanism || {};
const validation = DATA.validation || {};

const ACTION_LABEL = {
  define_discordance:"定义不一致", discover_subtypes:"亚型发现", select_slice:"切片复现",
  test_confounder:"混杂归因", sensitivity_analysis:"敏感性", profile_mechanism:"机制定位"
};

function keyOf(e){
  const m = (e.feedback && e.feedback.metrics) || {};
  const a = e.action;
  try{
    if(a==="define_discordance"){
      const s=m.discordant_group_size||{};
      return "不一致率 "+(s.discordance_rate*100).toFixed(1)+"% · 四组 "+JSON.stringify(s.by_group);
    }
    if(a==="select_slice"){
      const c=(m.cross_slice_reproducibility||{}).slice_intermediate_check||{};
      return "中间态 "+(c.holds?"成立":"不成立");
    }
    if(a==="test_confounder"){
      const d=(m.confounder_adjusted_change||{}).delta;
      return "中间态系数 Δ="+(d==null?"—":Number(d).toFixed(4));
    }
    if(a==="discover_subtypes"){
      const sil=(m.discordant_group_size||{}).silhouette;
      const p=(m.trajectory_separation||{}).mw_p;
      return "轮廓 "+(sil==null?"—":Number(sil).toFixed(2))+" · 轨迹 p="+(p==null?"—":Number(p).toFixed(3));
    }
    if(a==="sensitivity_analysis"){
      const c=(m.cross_slice_reproducibility||{}).coeff;
      return "复现系数 "+(c==null?"—":Number(c).toFixed(2));
    }
  }catch(err){}
  return "";
}

function badgeClass(e){
  if(e.action==="define_discordance" && e.params && e.params.missing_policy==="exclude") return "fix";
  const p=(e.feedback&&e.feedback.metrics&&e.feedback.metrics.trajectory_separation)||{};
  if(p.mw_p!=null && p.mw_p>=0.05) return "neg";
  return "";
}

function renderKpis(){
  const r2 = logs.find(e=>e.round_id===2);
  const s = r2 && r2.feedback.metrics.discordant_group_size;
  const kpis = [];
  if(s) kpis.push(["修正不一致率", (s.discordance_rate*100).toFixed(1)+"%"]);
  if(s) kpis.push(["不一致人数", s.n_discordant]);
  if(s) kpis.push(["PET+/Plasma−", s.by_group["PET+/Plasma−"]]);
  kpis.push(["探索轮次", logs.length]);
  kpis.push(["参照系", "4 层"]);
  document.getElementById("kpis").innerHTML = kpis.map(k=>
    '<div class="kpi"><b>'+k[1]+'</b><span>'+k[0]+'</span></div>').join("");
}

function renderTimeline(){
  const el = document.getElementById("timeline");
  el.innerHTML = logs.map(e=>{
    const bc = badgeClass(e);
    const node = bc==="fix" ? "node fix" : (bc==="neg" ? "node neg" : "node");
    const params = e.params ? JSON.stringify(e.params) : "{}";
    const ok = e.success===false ? "✗" : "✓";
    return '<div class="round"><span class="'+node+'"></span>'
      +'<div class="head"><span class="badge '+(bc==="fix"?"fix":bc==="neg"?"neg":"")+'">'
      +'R'+e.round_id+' · '+(ACTION_LABEL[e.action]||e.action)+'</span>'
      +'<span class="tag">state_update: '+e.state_update+'</span>'
      +'<span class="tag">'+ok+'</span>'
      +'<span class="tag">'+params+'</span></div>'
      +'<div class="rationale">'+e.rationale+'</div>'
      +'<div class="key">'+keyOf(e)+'</div>'
      +'<details><summary>完整 5 项反馈</summary><pre>'+JSON.stringify(e.feedback,null,1)+'</pre></details>'
      +'</div>';
  }).join("");
}

function fmt(v){
  if(v==null) return "—";
  if(typeof v==="number") return Number(v).toFixed(3);
  return String(v);
}

function renderRef(){
  const t = ref.trivial||{};
  const r = ref.random||{};
  const n = ref.no_intervention||{};
  const nt = ref.nontrivial||{};
  const cards = [
    {h:"① 平凡解", v:"仅PET d="+fmt(t.pet_only_binary&&t.pet_only_binary.cohens_d&&t.pet_only_binary.cohens_d.d),
     p:"单模态-仅 PET 分离 "+(t.pet_only_binary&&t.pet_only_binary.mw_p?fmt(t.pet_only_binary.mw_p):"—")+"；仅血浆 d="+fmt(t.plasma_only_binary&&t.plasma_only_binary.cohens_d&&t.plasma_only_binary.cohens_d.d)},
    {h:"② 随机参照", v:r.verdict||"—",
     p:"四组 H_obs="+fmt(r.H_obs)+" vs 置换 p95="+fmt(r.perm_dist&&r.perm_dist.p95)+" · perm_p="+fmt(r.perm_p)},
    {h:"③ 无干预参照", v:"中间态 p="+fmt(n.intermediate_vs_neg&&n.intermediate_vs_neg.p),
     p:"四组中位 "+JSON.stringify(n.four_group_medians)},
    {h:"④ Pyun baseline", v:"跨切片复现 "+fmt(nt.add_value_layer_2_cross_slice&&nt.add_value_layer_2_cross_slice.repro_coeff),
     p:"混杂校正 ΔG1="+fmt(nt.add_value_layer_3_confounder&&nt.add_value_layer_3_confounder.delta)}
  ];
  document.getElementById("ref").innerHTML = cards.map(c=>
    '<div class="ref-card"><h3>'+c.h+'</h3><div class="v">'+c.v+'</div><p>'+c.p+'</p></div>').join("");
}

function renderMechanisms(){
  const profiles = mechanism.profiles || {};
  const cards = Object.keys(profiles).map(group=>{
    const p = profiles[group] || {};
    const ranking = p.candidate_mechanism_ranking || [];
    const top = ranking.slice(0,3).map(x=>x.mechanism+"="+fmt(x.score)).join(" · ");
    const injury = p.injury_evidence || {};
    const gfap = injury.GFAP && injury.GFAP.effect ? injury.GFAP.effect.median_difference : null;
    const nfl = injury.NfL && injury.NfL.effect ? injury.NfL.effect.median_difference : null;
    const tau = p.tau_evidence && p.tau_evidence.effect ? p.tau_evidence.effect.median_difference : null;
    return '<div class="ref-card"><h3>'+group+'</h3><div class="v">'+top+'</div>'
      +'<p>tau Δ='+fmt(tau)+' · GFAP Δ='+fmt(gfap)+' · NfL Δ='+fmt(nfl)+'</p>'
      +'<p>对齐窗口 '+(p.alignment_window||"—")+' · n='+fmt(p.sample&&p.sample.n_target)+'</p></div>';
  });
  document.getElementById("mechanisms").innerHTML = cards.join("") || '<div class="card">暂无机制报告</div>';
  const split = validation.split || {};
  const ranks = validation.ranking_reproduction || {};
  document.getElementById("validation").innerHTML = '<p><b>Site-held-out validation：</b>训练 '+fmt(split.n_train)+'，验证 '+fmt(split.n_validation)+'，留出 site '+((split.holdout_sites||[]).join(", ")||"—")+'</p>'
    +Object.keys(ranks).map(group=>'<p><b>'+group+'：</b>Top-2 overlap '+fmt(ranks[group].top2_overlap)+'；训练 '+(ranks[group].train_top2||[]).join(" / ")+'；验证 '+(ranks[group].validation_top2||[]).join(" / ")+'</p>').join("");
}

function renderFindings(){
  document.getElementById("findings").innerHTML =
    '<p><b>轨迹：</b>PET−/Plasma+ 温和中间态（ΔADAS13/年 +0.46，vs 双阴 p=0.033）；PET+/Plasma− 与双阴无差异（+0.39 vs +0.21，p=0.94）——「影像先行型接近双阳」不成立。</p>'
    +'<p><b>合并症：</b>负担最高的是 PET−/Plasma+（心血管 34%/内分泌 25%），而非 PET+/Plasma−（修正后 10.8%/12.5%）。</p>'
    +'<p><b>tau：</b>两方向趋同（PET+/Plasma− 1.195 ≈ PET−/Plasma+ 1.191 ≈ 双阴 1.169，远低于双阳 1.387），均无显著 tau 启动。</p>'
    +'<div class="callout"><b>统一叙事：</b>排除缺失口径后，轨迹/合并症/tau 三条证据链收敛——两个不一致方向均接近「无 tau 启动」；PET−/Plasma+ 才是合并症驱动、值得继续追踪的组；PET+/Plasma− 无独立信号（负结果同样有价值）。</div>';
}

renderKpis();
renderTimeline();
renderRef();
renderFindings();
renderMechanisms();
</script>
</body>
</html>
"""


def main():
    logs = load_log()
    ref = load_ref()
    mechanism = load_json("mechanism_report.json")
    validation = load_json("validation_report.json")
    payload = json.dumps({"logs": logs, "ref": ref, "mechanism": mechanism, "validation": validation}, ensure_ascii=False)
    html = HTML.replace("__DATA__", payload)
    out = ROOT / "dashboard.html"
    out.write_text(html, encoding="utf-8")
    print(f"仪表盘已生成：{out}（{out.stat().st_size // 1024} KB）—— 双击用浏览器打开即可")


if __name__ == "__main__":
    main()
