# -*- coding: utf-8 -*-
"""webapp/server.py — GOAI 本地探索台（零第三方依赖：标准库后端 + 原生前端）

功能：
  - 总览：四组结构 / 认知轨迹 / tau 负担实时统计（pandas 现算 + 60s 缓存）
  - LLM 探索：OpenAI 兼容接口驱动 ExplorationAgent 闭环（无 key 自动降级 MockPolicy）
  - 数据替换：上传任意管线表 CSV 覆盖 data/raw/merged/，一键重跑管线
  - 报告：浏览 data/processed/*.txt 聚合报告

运行：python webapp/server.py  →  http://127.0.0.1:8765
"""
import json
import os
import re
import subprocess
import sys
import threading
import time
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs

ROOT = Path(__file__).resolve().parents[1]
WEBAPP_DIR = ROOT / "webapp"
sys.path.insert(0, str(ROOT / "scripts"))

import config as C  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from trajectory import per_patient_change  # noqa: E402

HOST, PORT = "127.0.0.1", 8765
TABLE_FILES = {Path(f).name for f in C.FILES.values()}  # 可替换的表文件名白名单（basename）
GROUP_LABELS = {0: "PET−/Plasma−", 1: "PET−/Plasma+", 2: "PET+/Plasma−", 3: "PET+/Plasma+"}

# ---------------- 后台任务 ----------------
# 注意：start_task 在持锁状态下调用 _set_task（内部也拿锁），必须用可重入锁，否则死锁。
_task_lock = threading.RLock()
_task = {"status": "idle", "kind": None, "output": "", "started": None, "finished": None}


def _set_task(status=None, kind=None, output=None, append=None):
    with _task_lock:
        if status is not None:
            _task["status"] = status
        if kind is not None:
            _task["kind"] = kind
        if append is not None:
            _task["output"] += append
        if output is not None:
            _task["output"] = output
        if status == "running":
            _task["started"] = datetime.now().isoformat(timespec="seconds")
        if status in ("done", "error"):
            _task["finished"] = datetime.now().isoformat(timespec="seconds")


def _get_task():
    with _task_lock:
        return dict(_task)


# ---------------- 统计（缓存 60s，任务完成/上传后失效） ----------------
_stats_lock = threading.Lock()
_stats_cache = {"at": 0.0, "data": None}


def invalidate_stats():
    with _stats_lock:
        _stats_cache["at"] = 0.0


def _wide() -> pd.DataFrame:
    df = pd.read_csv(C.PROC_DIR / "subjects_wide.csv")
    return df


def _groups(df: pd.DataFrame) -> pd.DataFrame:
    sub = df.dropna(subset=["PET_STATUS", "PLASMA_STATUS"]).copy()
    sub["GROUP"] = sub["PET_STATUS"].astype(int) * 2 + sub["PLASMA_STATUS"].astype(int)
    return sub


def _tau_medians(sub: pd.DataFrame) -> dict:
    tau = pd.read_csv(C.RAW_DIR / C.FILES["tau"], low_memory=False)
    tau = tau[(tau["TRACER"] == "FTP") & (tau["qc_flag"].isin([1, 2]))].copy()
    tau["SCANDATE"] = pd.to_datetime(tau["SCANDATE"], errors="coerce")
    tau = tau.dropna(subset=["SCANDATE", C.FIELDS["tau_suvr"]])
    tau = tau.sort_values(["RID", "SCANDATE"]).groupby("RID", as_index=False).first()
    m = sub.merge(tau[["RID", C.FIELDS["tau_suvr"]]], on="RID", how="left")
    out = {}
    for g in sorted(sub.GROUP.unique()):
        v = m.loc[m.GROUP == g, C.FIELDS["tau_suvr"]].dropna()
        out[g] = {"n": int(len(v)), "median": float(v.median()) if len(v) else None,
                  "q1": float(v.quantile(.25)) if len(v) else None,
                  "q3": float(v.quantile(.75)) if len(v) else None,
                  "min": float(v.min()) if len(v) else None,
                  "max": float(v.max()) if len(v) else None}
    return out


def _trajectory_medians(sub: pd.DataFrame) -> dict:
    cdr = pd.read_csv(C.RAW_DIR / C.FILES["cdr"], low_memory=False)
    adas = pd.read_csv(C.RAW_DIR / C.FILES["adas"], low_memory=False)
    c = per_patient_change(cdr, C.FIELDS["date"]["cdr"], C.FIELDS["cdrsb"], "D_CDRSB", "YRS_CDR")
    a = per_patient_change(adas, C.FIELDS["date"]["adas"], C.FIELDS["adas13"], "D_ADAS13", "YRS_ADAS")
    m = sub.merge(c, on="RID", how="left").merge(a, on="RID", how="left")
    m["D_CDRSB_yr"] = m["D_CDRSB"] / m["YRS_CDR"]
    m["D_ADAS13_yr"] = m["D_ADAS13"] / m["YRS_ADAS"]
    out = {}
    for g in sorted(sub.GROUP.unique()):
        row = {}
        for col in ("D_ADAS13_yr", "D_CDRSB_yr"):
            v = m.loc[m.GROUP == g, col].dropna()
            row[col] = {"n": int(len(v)), "median": float(v.median()) if len(v) else None}
        out[g] = row
    return out


def _comorbidity_rates(sub: pd.DataFrame) -> dict:
    mh = pd.read_csv(C.RAW_DIR / C.FILES["medhist"], low_memory=False)
    mh = mh.sort_values(["RID", "VISDATE"]).groupby("RID", as_index=False).first()
    m = sub.merge(mh[["RID", "MH4CARD", "MH9ENDO"]], on="RID", how="left")
    out = {}
    for g in sorted(sub.GROUP.unique()):
        d = m[m.GROUP == g]
        out[g] = {"cardio": float(C.normalize_yes_no(d["MH4CARD"]).fillna(0).mean()),
                  "endo": float(C.normalize_yes_no(d["MH9ENDO"]).fillna(0).mean())}
    return out


def compute_stats() -> dict:
    with _stats_lock:
        if _stats_cache["data"] is not None and time.time() - _stats_cache["at"] < 60:
            return _stats_cache["data"]
    df = _wide()
    sub = _groups(df)
    n_total = int(len(df))
    traj, tau, com = _trajectory_medians(sub), _tau_medians(sub), _comorbidity_rates(sub)
    groups = []
    for g in sorted(sub.GROUP.unique()):
        d = sub[sub.GROUP == g]
        groups.append({
            "group": int(g), "label": GROUP_LABELS[g], "n": int(len(d)),
            "age_median": float(d.AGE.median()) if d.AGE.notna().any() else None,
            "adas13_yr": traj[g]["D_ADAS13_yr"], "cdrsb_yr": traj[g]["D_CDRSB_yr"],
            "tau": tau[g], "cardio": com[g]["cardio"], "endo": com[g]["endo"],
        })
    n_disc = sum(x["n"] for x in groups if x["group"] in (1, 2))
    data = {"n_total": n_total, "groups": groups,
            "discordance_rate": n_disc / max(1, sum(x["n"] for x in groups)),
            "generated_at": datetime.now().isoformat(timespec="seconds")}
    with _stats_lock:
        _stats_cache["at"], _stats_cache["data"] = time.time(), data
    return data


# ---------------- 文件与报告 ----------------
def tables_info() -> list:
    out = []
    for fname in sorted(TABLE_FILES):
        p = C.RAW_DIR / fname
        if not p.exists():
            out.append({"file": fname, "exists": False, "rows": None, "mtime": None})
            continue
        try:
            rows = int(sum(1 for _ in open(p, "rb"))) - 1
        except OSError:
            rows = None
        out.append({"file": fname, "exists": True, "rows": rows,
                    "mtime": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")})
    return out


def reports_list() -> list:
    return sorted(p.name for p in (C.PROC_DIR).glob("*.txt"))


def llm_info() -> dict:
    p = ROOT / "llm_config.json"
    cfg = {}
    if p.exists():
        try:
            cfg = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            cfg = {}
    return {"configured": bool(cfg.get("api_key")), "model": cfg.get("model", ""),
            "base_url": cfg.get("base_url", "")}


def explore_summary() -> dict:
    p = C.PROC_DIR / "exploration_log.jsonl"
    if not p.exists():
        return {"rounds": 0, "latest": []}
    logs = []
    for line in p.read_text(encoding="utf-8").splitlines():
        try:
            logs.append(json.loads(line))
        except Exception:
            continue
    latest = []
    for e in logs[-6:]:
        fb = e.get("feedback", {})
        latest.append({"round": e.get("round_id"), "action": e.get("action"),
                       "params": e.get("params", {}), "rationale": e.get("rationale", ""),
                       "state_update": e.get("state_update", ""),
                       "ok": "error" not in fb,
                       "note": fb.get("interpretation", fb.get("error", ""))[:200]})
    return {"rounds": len(logs), "latest": latest}


# ---------------- multipart 解析（标准库无现成） ----------------
def parse_multipart(body: bytes, content_type: str) -> list:
    m = re.search(r"boundary=([^;]+)", content_type)
    if not m:
        raise ValueError("无 boundary")
    boundary = m.group(1).strip().strip('"').encode()
    parts, marker = [], b"--" + boundary
    for chunk in body.split(marker)[1:]:
        if chunk.startswith(b"--"):
            break
        head, _, content = chunk.partition(b"\r\n\r\n")
        head_txt = head.decode("utf-8", "replace")
        name_m = re.search(r'name="([^"]+)"', head_txt)
        file_m = re.search(r'filename="([^"]*)"', head_txt)
        if name_m:
            parts.append({"name": name_m.group(1),
                          "filename": file_m.group(1) if file_m else None,
                          "content": content[:-2] if content.endswith(b"\r\n") else content})
    return parts


# ---------------- 任务执行 ----------------
PIPELINE_STEPS = ["preprocess.py", "calibrate.py", "define_discordance.py",
                  "trajectory.py", "slice_analysis.py", "cluster_pm.py",
                  "tau_analysis.py", "tau_analysis_corrected.py", "corrected_reports.py"]
PIPELINE_SHORT = ["preprocess.py", "calibrate.py", "define_discordance.py"]


def _run_scripts(scripts: list, label: str):
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    for s in scripts:
        _set_task(append=f"\n$ python scripts/{s}\n")
        try:
            r = subprocess.run([sys.executable, str(ROOT / "scripts" / s)],
                               capture_output=True, text=True, cwd=ROOT, env=env, timeout=1800)
        except subprocess.TimeoutExpired:
            _set_task(status="error", append="\n[超时] 该步骤超过 30 分钟，中止。\n")
            return
        out = (r.stdout or "") + (r.stderr or "")
        _set_task(append=out[-3000:])
        if r.returncode != 0:
            _set_task(status="error", append=f"\n[失败] {s} 退出码 {r.returncode}\n")
            return
    _set_task(append=f"\n[完成] {label}\n")
    invalidate_stats()
    _set_task(status="done")


def _run_explore(rounds: int, api_key: str, base_url: str, model: str):
    from agent import ExplorationAgent, OpenAICompatLLM
    from environment import ExplorationEnvironment
    _set_task(append=f"\n[探索] 启动 {rounds} 轮闭环（LLM={'有' if api_key else 'MockPolicy'}）\n")
    env = ExplorationEnvironment()
    llm = OpenAICompatLLM(api_key=api_key or "", base_url=base_url or None, model=model or None) \
        if api_key else None
    agent = ExplorationAgent(env, llm=llm)
    logs, path = agent.run(n_rounds=rounds, verbose=False)
    _set_task(append=f"[探索] 日志已写 {path}（{len(logs)} 轮）\n")
    invalidate_stats()
    _set_task(status="done")


def start_task(kind: str, **kw):
    with _task_lock:
        if _task["status"] == "running":
            return False
        _set_task(status="running", kind=kind, output=f"任务启动：{kind}\n")
    if kind == "pipeline":
        threading.Thread(target=_run_scripts, args=(PIPELINE_STEPS, "完整管线"), daemon=True).start()
    elif kind == "quick":
        threading.Thread(target=_run_scripts, args=(PIPELINE_SHORT, "快速管线（preprocess→四组）"), daemon=True).start()
    elif kind == "explore":
        threading.Thread(target=_run_explore, args=(kw["rounds"], kw.get("api_key", ""),
                                                    kw.get("base_url", ""), kw.get("model", "")), daemon=True).start()
    return True


# ---------------- HTTP ----------------
class Handler(BaseHTTPRequestHandler):
    server_version = "GOAIExplorer/1.0"

    def log_message(self, fmt, *args):
        pass

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, ensure_ascii=False).encode("utf-8")
        elif isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        if u.path in ("/", "/index.html"):
            return self._send(200, (WEBAPP_DIR / "index.html").read_text(encoding="utf-8"), "text/html; charset=utf-8")
        if u.path == "/style.css":
            return self._send(200, (WEBAPP_DIR / "style.css").read_text(encoding="utf-8"), "text/css; charset=utf-8")
        if u.path == "/app.js":
            return self._send(200, (WEBAPP_DIR / "app.js").read_text(encoding="utf-8"), "application/javascript; charset=utf-8")
        if u.path == "/api/state":
            return self._send(200, {"stats": compute_stats(), "tables": tables_info(),
                                    "reports": reports_list(), "llm": llm_info(),
                                    "task": _get_task(), "explore": explore_summary()})
        if u.path == "/api/report":
            name = parse_qs(u.query).get("name", [""])[0]
            if not re.fullmatch(r"[\w\-]+\.txt", name):
                return self._send(400, {"error": "非法报告名"})
            p = C.PROC_DIR / name
            if not p.exists():
                return self._send(404, {"error": "报告不存在"})
            return self._send(200, p.read_text(encoding="utf-8"), "text/plain; charset=utf-8")
        if u.path == "/api/task":
            return self._send(200, _get_task())
        return self._send(404, {"error": "not found"})

    def do_POST(self):
        u = urlparse(self.path)
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length) if length else b""

        if u.path == "/api/upload":
            try:
                parts = parse_multipart(body, self.headers.get("Content-Type", ""))
            except ValueError:
                return self._send(400, {"error": "需要 multipart/form-data"})
            saved, skipped = [], []
            for part in parts:
                fname = part.get("filename", "")
                if fname not in TABLE_FILES:
                    skipped.append(f"{fname}（不在白名单）")
                    continue
                if len(part["content"]) > 300 * 1024 * 1024:
                    skipped.append(f"{fname}（超过 300MB 上限）")
                    continue
                (C.RAW_DIR / fname).write_bytes(part["content"])
                saved.append(fname)
            invalidate_stats()
            return self._send(200, {"saved": saved, "skipped": skipped,
                                    "hint": "替换后请到「数据」页重跑管线（快速或完整）"})

        if u.path == "/api/llm_config":
            try:
                cfg = json.loads(body.decode("utf-8"))
            except Exception:
                return self._send(400, {"error": "JSON 解析失败"})
            if cfg.get("remember"):
                (ROOT / "llm_config.json").write_text(
                    json.dumps({"api_key": cfg.get("api_key", ""), "base_url": cfg.get("base_url", ""),
                                "model": cfg.get("model", "")}, ensure_ascii=False, indent=2),
                    encoding="utf-8")
            return self._send(200, {"ok": True})

        if u.path == "/api/task":
            try:
                req = json.loads(body.decode("utf-8"))
            except Exception:
                return self._send(400, {"error": "JSON 解析失败"})
            kind = req.get("kind")
            if kind == "explore":
                rounds = max(1, min(int(req.get("rounds", 8)), 20))
                ok = start_task("explore", rounds=rounds, api_key=req.get("api_key", ""),
                                base_url=req.get("base_url", ""), model=req.get("model", ""))
            elif kind in ("pipeline", "quick"):
                ok = start_task(kind)
            else:
                return self._send(400, {"error": f"未知任务类型: {kind}"})
            if not ok:
                return self._send(409, {"error": "已有任务在运行"})
            return self._send(200, {"ok": True})

        return self._send(404, {"error": "not found"})


def main():
    httpd = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f"GOAI 探索台已启动：http://{HOST}:{PORT}  （Ctrl+C 退出）")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已退出")


if __name__ == "__main__":
    main()
