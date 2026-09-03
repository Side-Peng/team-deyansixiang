/* GOAI 探索台前端逻辑（零依赖：原生 JS + Canvas） */
"use strict";

const GROUP_COLORS = { 0: "#64748b", 1: "#fbbf24", 2: "#7aa2ff", 3: "#f472b6" };
const $ = (s) => document.querySelector(s);

/* ---------------- 基础 ---------------- */
async function api(path, opts) {
  const r = await fetch(path, opts);
  if (!r.ok) {
    let msg = r.statusText;
    try { msg = (await r.json()).error || msg; } catch (e) {}
    throw new Error(msg);
  }
  return r.json();
}

function fmtPct(x) { return x == null ? "—" : (x * 100).toFixed(1) + "%"; }
function fmtNum(x, d = 2) { return x == null ? "—" : Number(x).toFixed(d); }

let state = null;
let pollTimer = null;

/* ---------------- 总览渲染 ---------------- */
function renderKPIs() {
  const s = state.stats;
  const g1 = s.groups.find(g => g.group === 1);
  const g2 = s.groups.find(g => g.group === 2);
  const g3 = s.groups.find(g => g.group === 3);
  const kpis = [
    { label: "不一致率（修正口径）", value: fmtPct(s.discordance_rate), sub: `${s.groups.map(g => g.n).join(" / ")}（四组）` },
    { label: "基线总人数", value: s.n_total, sub: "subjects_wide.csv" },
    { label: "中间态 P−/P+ 轨迹", value: g1 ? "+" + fmtNum(g1.adas13_yr.median) : "—", sub: `ΔADAS13/年（n=${g1?.adas13_yr.n ?? 0}）` },
    { label: "影像先行 PET+/P− 轨迹", value: g2 ? "+" + fmtNum(g2.adas13_yr.median) : "—", sub: `ΔADAS13/年（n=${g2?.adas13_yr.n ?? 0}）` },
    { label: "tau SUVR：P−/P+ vs PET+/P−", value: g1 && g2 ? `${fmtNum(g1.tau.median)} / ${fmtNum(g2.tau.median)}` : "—", sub: `双阳 ${g3 ? fmtNum(g3.tau.median) : "—"}（两方向趋同？）` },
    { label: "探索闭环轮数", value: state.explore.rounds, sub: state.explore.rounds ? "exploration_log.jsonl" : "尚未运行" },
  ];
  $("#kpis").innerHTML = kpis.map(k =>
    `<div class="kpi"><div class="k-label">${k.label}</div><div class="k-value">${k.value}</div><div class="k-sub">${k.sub}</div></div>`
  ).join("");
}

function drawChart(canvasId, renderFn) {
  const cv = $(canvasId);
  const dpr = window.devicePixelRatio || 1;
  const w = cv.parentElement.clientWidth - 40;
  const h = +cv.getAttribute("height");
  cv.width = w * dpr; cv.height = h * dpr;
  cv.style.width = w + "px";
  const ctx = cv.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);
  renderFn(ctx, w, h);
}

/* 轨迹：分组条形图（4 组 × 2 结局） */
function renderTrajectory() {
  const s = state.stats;
  drawChart("#chartTrajectory", (ctx, w, h) => {
    const padL = 52, padB = 44, padT = 16;
    const plotW = w - padL - 16, plotH = h - padT - padB;
    const outcomes = [["adas13_yr", "#7aa2ff", "ΔADAS13/年"], ["cdrsb_yr", "#5eead4", "ΔCDRSB/年"]];
    const vals = outcomes.flatMap(([k]) => s.groups.map(g => g[k].median ?? 0));
    const maxV = Math.max(0.4, ...vals) * 1.15;
    const minV = Math.min(0, ...vals) * 1.1;
    const y = v => padT + plotH - (v - minV) / (maxV - minV) * plotH;

    ctx.strokeStyle = "rgba(148,180,255,.15)"; ctx.fillStyle = "#8b95ab"; ctx.font = "11px Consolas";
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (let t = 0; t <= 4; t++) {
      const v = minV + (maxV - minV) * t / 4, yy = y(v);
      ctx.beginPath(); ctx.moveTo(padL, yy); ctx.lineTo(w - 16, yy); ctx.stroke();
      ctx.fillText(v.toFixed(1), padL - 8, yy);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    const band = plotW / s.groups.length;
    s.groups.forEach((g, gi) => {
      const x0 = padL + gi * band;
      ctx.fillStyle = "#c7d2e8"; ctx.font = "12px Segoe UI";
      ctx.fillText(g.label.replace("/", " /\n").split("\n").join(" / "), x0 + band / 2, h - padB + 8);
      ctx.fillStyle = GROUP_COLORS[g.group];
      const bw = band * 0.30;
      outcomes.forEach(([k, color], oi) => {
        const v = g[k].median ?? 0;
        const bx = x0 + band / 2 - bw + oi * bw * 0.75;
        ctx.fillStyle = color;
        rounded(ctx, bx, y(Math.max(0, v)), bw * 0.7, Math.abs(y(v) - y(0)));
        ctx.fillStyle = color; ctx.font = "10px Consolas"; ctx.textAlign = "center";
        ctx.fillText(v.toFixed(2), bx + bw * 0.35, y(Math.max(0, v)) - 14);
      });
    });
    ctx.fillStyle = "#8b95ab"; ctx.textAlign = "left"; ctx.font = "11px Segoe UI";
    outcomes.forEach(([k, color, name], oi) => {
      ctx.fillStyle = color; ctx.fillRect(16, 6 + oi * 16, 9, 9);
      ctx.fillStyle = "#8b95ab"; ctx.fillText(name, 30, 11 + oi * 16);
    });
  });
}

/* tau：四组箱线图 */
function renderTau() {
  const s = state.stats;
  drawChart("#chartTau", (ctx, w, h) => {
    const padL = 52, padB = 40, padT = 20;
    const plotW = w - padL - 16, plotH = h - padT - padB;
    const vals = s.groups.flatMap(g => g.tau.n ? [g.tau.min, g.tau.q1, g.tau.median, g.tau.q3, g.tau.max] : []);
    if (!vals.length) { ctx.fillStyle = "#8b95ab"; ctx.fillText("无 tau 数据", w / 2, h / 2); return; }
    const lo = Math.min(...vals) * 0.96, hi = Math.max(...vals) * 1.04;
    const y = v => padT + plotH - (v - lo) / (hi - lo) * plotH;

    ctx.strokeStyle = "rgba(148,180,255,.15)"; ctx.fillStyle = "#8b95ab"; ctx.font = "11px Consolas";
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (let t = 0; t <= 4; t++) {
      const v = lo + (hi - lo) * t / 4, yy = y(v);
      ctx.beginPath(); ctx.moveTo(padL, yy); ctx.lineTo(w - 16, yy); ctx.stroke();
      ctx.fillText(v.toFixed(2), padL - 8, yy);
    }
    // 双阳参考线
    const g3 = s.groups.find(g => g.group === 3);
    if (g3?.tau.median != null) {
      ctx.strokeStyle = "rgba(244,114,182,.55)"; ctx.setLineDash([5, 5]);
      ctx.beginPath(); ctx.moveTo(padL, y(g3.tau.median)); ctx.lineTo(w - 16, y(g3.tau.median)); ctx.stroke();
      ctx.setLineDash([]); ctx.fillStyle = "#f472b6"; ctx.textAlign = "left";
      ctx.fillText(`双阳 ${g3.tau.median.toFixed(2)}`, w - 96, y(g3.tau.median) - 7);
    }
    const band = plotW / s.groups.length;
    ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.font = "12px Segoe UI";
    s.groups.forEach((g, gi) => {
      const x0 = padL + gi * band, cx = x0 + band / 2;
      ctx.fillStyle = "#c7d2e8";
      ctx.fillText(g.label, cx, h - padB + 8);
      ctx.fillStyle = GROUP_COLORS[g.group];
      ctx.fillText("n=" + (g.tau.n ?? 0), cx, h - padB + 24);
      if (!g.tau.n) return;
      const t = g.tau, bw = band * 0.42;
      ctx.strokeStyle = GROUP_COLORS[g.group]; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(cx - bw / 2, y(t.min)); ctx.lineTo(cx + bw / 2, y(t.min)); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(cx, y(t.min)); ctx.lineTo(cx, y(t.q1)); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(cx, y(t.q3)); ctx.lineTo(cx, y(t.max)); ctx.stroke();
      ctx.fillStyle = "rgba(122,162,255,.16)"; ctx.strokeStyle = GROUP_COLORS[g.group]; ctx.lineWidth = 1.5;
      ctx.fillRect(cx - bw / 2, y(t.q3), bw, y(t.q1) - y(t.q3));
      ctx.strokeRect(cx - bw / 2, y(t.q3), bw, y(t.q1) - y(t.q3));
      ctx.fillStyle = GROUP_COLORS[g.group]; ctx.font = "10px Consolas";
      ctx.fillText(t.median.toFixed(2), cx, y(t.median) - 12);
    });
  });
}

function renderTimeline() {
  const items = state.explore.latest;
  if (!items.length) { $("#timeline").innerHTML = '<p class="empty">暂无探索记录 —— 到「LLM 探索」页启动闭环</p>'; return; }
  $("#timeline").innerHTML = items.map(e => `
    <div class="tl-item">
      <div class="tl-round">R${e.round}</div>
      <div class="tl-action">${e.action}</div>
      <div class="tl-rationale" title="${esc(e.rationale)}">${esc(e.rationale)}</div>
      <span class="badge ${e.ok ? "ok" : "fail"}">${e.ok ? "✓" : "✗"} ${esc(e.state_update || "")}</span>
    </div>`).join("");
}

function rounded(ctx, x, y, w, h, r = 3) {
  ctx.beginPath();
  ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r);
  ctx.arcTo(x + w, y + h, x, y + h, r); ctx.arcTo(x, y + h, x, y, r);
  ctx.arcTo(x, y, x + w, y, r); ctx.closePath(); ctx.fill();
}

function esc(s) { return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])); }

/* ---------------- 总览刷新 ---------------- */
async function refreshState(force) {
  try {
    state = await api("/api/state");
    renderKPIs(); renderTrajectory(); renderTau(); renderTimeline();
    renderTables(); renderReports(); renderLLM();
    renderTask();
    setStatus("ok", "已连接 · " + state.stats.generated_at.slice(11));
  } catch (e) {
    setStatus("err", "连接失败");
  }
  return state;
}

function setStatus(cls, text) {
  const b = $("#statusBadge");
  b.className = "status " + cls;
  b.innerHTML = `<span class="dot"></span>${text}`;
}

/* ---------------- LLM 面板 ---------------- */
function renderLLM() {
  const llm = state.llm;
  $("#llmBase").value = $("#llmBase").value || llm.base_url || "";
  $("#llmModel").value = $("#llmModel").value || llm.model || "";
}

$("#llmForm").addEventListener("submit", async (e) => {
  e.preventDefault();
  const body = {
    api_key: $("#llmKey").value.trim(),
    base_url: $("#llmBase").value.trim(),
    model: $("#llmModel").value.trim(),
    remember: $("#llmRemember").checked,
  };
  await api("/api/llm_config", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  const rounds = +$("#llmRounds").value;
  try {
    await api("/api/task", { method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ kind: "explore", rounds, ...body }) });
  } catch (err) {
    alert("启动失败：" + err.message);
    return;
  }
  $("#exploreBtn").disabled = true;
  $("#taskState").textContent = "running";
  pollTask();
});

/* ---------------- 任务轮询 ---------------- */
function renderTask() {
  const t = state.task;
  const out = $("#taskOutput");
  if (t.status === "idle" && !t.output) { out.textContent = "—"; return; }
  out.textContent = t.output || "—";
  out.scrollTop = out.scrollHeight;
  $("#taskState").textContent = t.status + (t.kind ? " · " + t.kind : "");
  const busy = t.status === "running";
  $("#exploreBtn").disabled = busy;
  $("#btnPipeline").disabled = busy;
  $("#btnQuick").disabled = busy;
  setStatus(busy ? "busy" : t.status === "error" ? "err" : "ok",
    busy ? "任务运行中 · " + (t.kind || "") : "已连接");
}

function pollTask() {
  clearInterval(pollTimer);
  pollTimer = setInterval(async () => {
    try {
      const t = await api("/api/task");
      const out = $("#taskOutput");
      out.textContent = t.output || "—";
      out.scrollTop = out.scrollHeight;
      $("#taskState").textContent = t.status + (t.kind ? " · " + t.kind : "");
      if (t.status !== "running") {
        clearInterval(pollTimer);
        const busy = false;
        $("#exploreBtn").disabled = busy; $("#btnPipeline").disabled = busy; $("#btnQuick").disabled = busy;
        $("#taskOutput").hidden = false;
        await refreshState();
        setStatus(t.status === "error" ? "err" : "ok", t.status === "error" ? "任务失败" : "任务完成");
        $("#taskOutput").textContent = t.output;
      }
    } catch (e) { clearInterval(pollTimer); }
  }, 1200);
}

/* ---------------- 数据面板 ---------------- */
function renderTables() {
  const rows = state.tables.map(t => `
    <tr><td>${t.exists ? esc(t.file) : "<s style='color:#f87171'>" + esc(t.file) + "</s>"}</td>
    <td>${t.rows == null ? "缺失" : t.rows.toLocaleString()}</td>
    <td>${t.mtime ? t.mtime.replace("T", " ") : "—"}</td></tr>`).join("");
  $("#tablesTable tbody").innerHTML = rows;
}

async function startPipeline(kind) {
  try {
    await api("/api/task", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ kind }) });
  } catch (e) { alert("启动失败：" + e.message); return; }
  $("#pipelineOutput").hidden = false;
  $("#pipelineOutput").textContent = "任务启动…";
  $("#btnPipeline").disabled = true; $("#btnQuick").disabled = true;
  pollTask();
}

$("#btnPipeline").addEventListener("click", () => startPipeline("pipeline"));
$("#btnQuick").addEventListener("click", () => startPipeline("quick"));

/* 上传 */
const dz = $("#dropZone"), fi = $("#fileInput");
dz.addEventListener("click", () => fi.click());
fi.addEventListener("change", () => uploadFiles(fi.files));
dz.addEventListener("dragover", e => { e.preventDefault(); dz.classList.add("over"); });
dz.addEventListener("dragleave", () => dz.classList.remove("over"));
dz.addEventListener("drop", e => { e.preventDefault(); dz.classList.remove("over"); uploadFiles(e.dataTransfer.files); });

async function uploadFiles(files) {
  if (!files.length) return;
  const fd = new FormData();
  for (const f of files) fd.append("files", f, f.name);
  const r = $("#uploadResult");
  r.textContent = "上传中…";
  try {
    const res = await api("/api/upload", { method: "POST", body: fd });
    r.className = "hint " + (res.skipped.length ? "warn" : "ok");
    r.textContent = `已替换：${res.saved.join("、") || "无"}${res.skipped.length ? "　跳过：" + res.skipped.join("；") : ""}`;
    if (res.saved.length) r.textContent += "　→ 请点击「快速重跑」或「完整重跑管线」让新数据生效。";
    refreshState();
  } catch (e) {
    r.className = "hint warn";
    r.textContent = "上传失败：" + e.message;
  }
}

/* ---------------- 报告面板 ---------------- */
function renderReports() {
  $("#reportList").innerHTML = state.reports.map(n =>
    `<li data-name="${esc(n)}">${esc(n)}</li>`).join("");
  $("#reportList").querySelectorAll("li").forEach(li =>
    li.addEventListener("click", async () => {
      $("#reportList").querySelectorAll("li").forEach(x => x.classList.remove("active"));
      li.classList.add("active");
      $("#reportTitle").textContent = li.dataset.name;
      $("#reportBody").textContent = await (await fetch(`/api/report?name=${encodeURIComponent(li.dataset.name)}`)).text();
    }));
}

/* ---------------- Tab 切换 ---------------- */
document.querySelectorAll(".tab").forEach(t =>
  t.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach(x => x.classList.remove("active"));
    document.querySelectorAll(".panel").forEach(x => x.classList.remove("active"));
    t.classList.add("active");
    $("#tab-" + t.dataset.tab).classList.add("active");
    if (t.dataset.tab === "overview") { renderTrajectory(); renderTau(); }
  }));

/* 窗口缩放时重绘图表 */
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if ($("#tab-overview").classList.contains("active")) { renderTrajectory(); renderTau(); }
  }, 200);
});

/* ---------------- 启动 ---------------- */
(async function init() {
  state = await refreshState();
  if (state.task.status === "running") pollTask();
})();
