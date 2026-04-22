#!/usr/bin/env python3
"""
bioBakery Nextflow Visualization Workflow — Monitoring Dashboard

A lightweight web dashboard for real-time monitoring of running biobakery
Nextflow visualization jobs. Shows pipeline progress, Nextflow trace stats,
SLURM queue state, and a live log tail.

Requirements:
    pip install fastapi uvicorn psutil

Usage:
    python monitoring/dashboard.py [--port 8050] [--output vis_output]
    open http://localhost:8050
"""

import argparse
import csv
import json
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

try:
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse, JSONResponse
    import uvicorn
    HAS_FASTAPI = True
except ImportError:
    HAS_FASTAPI = False


# ── Helpers ───────────────────────────────────────────────────────────────────

def read_trace(trace_file: Path) -> list[dict]:
    """Parse Nextflow trace.tsv into a list of task dicts."""
    if not trace_file.exists():
        return []
    rows = []
    with open(trace_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def read_log_tail(log_file: Path, n: int = 60) -> str:
    """Return the last n lines of .nextflow.log."""
    if not log_file.exists():
        return "Log file not found."
    try:
        lines = log_file.read_text(errors="replace").splitlines()
        return "\n".join(lines[-n:])
    except Exception as e:
        return f"Error reading log: {e}"


def pipeline_summary(trace_rows: list[dict]) -> dict:
    """Count tasks by status from trace rows."""
    counts = {"queued": 0, "running": 0, "completed": 0, "failed": 0, "cached": 0}
    for row in trace_rows:
        status = row.get("status", "").lower()
        if status == "completed":
            counts["completed"] += 1
        elif status == "failed":
            counts["failed"] += 1
        elif status in ("running", "submitted"):
            counts["running"] += 1
        elif status == "cached":
            counts["cached"] += 1
        else:
            counts["queued"] += 1
    counts["total"] = sum(counts.values())
    return counts


def slurm_jobs() -> list[dict]:
    """Run squeue and return current SLURM jobs as a list of dicts."""
    try:
        result = subprocess.run(
            ["squeue", "--format=%i,%j,%u,%T,%M,%l,%R,%C,%m",
             "--noheader", "--me"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode != 0:
            return []
        jobs = []
        headers = ["job_id","name","user","state","time","time_limit",
                   "reason","cpus","min_mem"]
        for line in result.stdout.strip().splitlines():
            parts = line.split(",", maxsplit=len(headers) - 1)
            jobs.append(dict(zip(headers, parts)))
        return jobs
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return []


def system_resources() -> dict:
    """Return current host CPU and memory usage."""
    if not HAS_PSUTIL:
        return {"cpu_pct": "N/A", "mem_pct": "N/A", "mem_used_gb": "N/A"}
    mem = psutil.virtual_memory()
    return {
        "cpu_pct": psutil.cpu_percent(interval=0.1),
        "mem_pct": mem.percent,
        "mem_used_gb": round(mem.used / 1e9, 1),
        "mem_total_gb": round(mem.total / 1e9, 1),
    }


# ── HTML template (single-file, no external CDN needed for air-gapped HPC) ───

DASHBOARD_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta http-equiv="refresh" content="15">
<title>bioBakery Vis — Nextflow Monitor</title>
<style>
  :root {{
    --bg: #1a1a2e; --panel: #16213e; --accent: #0f3460;
    --green: #4caf50; --red: #f44336; --yellow: #ff9800;
    --blue: #2196f3; --gray: #9e9e9e; --text: #e0e0e0;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: 'Courier New', monospace; font-size: 14px; }}
  header {{ background: var(--accent); padding: 16px 24px; display: flex; align-items: center; gap: 16px; }}
  header h1 {{ font-size: 20px; font-weight: bold; color: #fff; }}
  header .subtitle {{ color: #aaa; font-size: 12px; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; padding: 16px; }}
  .panel {{ background: var(--panel); border-radius: 8px; padding: 16px; border: 1px solid var(--accent); }}
  .panel h2 {{ font-size: 13px; color: #aaa; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 12px; }}
  .stat-row {{ display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px solid #223; }}
  .stat-row:last-child {{ border-bottom: none; }}
  .badge {{ padding: 2px 8px; border-radius: 10px; font-size: 11px; font-weight: bold; }}
  .badge-green {{ background: #1b5e20; color: var(--green); }}
  .badge-red {{ background: #b71c1c; color: #ff6666; }}
  .badge-yellow {{ background: #e65100; color: #ffd54f; }}
  .badge-blue {{ background: #0d47a1; color: #90caf9; }}
  .badge-gray {{ background: #333; color: var(--gray); }}
  .progress-bar {{ background: #333; border-radius: 4px; height: 10px; margin: 8px 0; }}
  .progress-fill {{ height: 10px; border-radius: 4px; background: var(--green); transition: width 0.3s; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  th {{ background: var(--accent); padding: 6px 8px; text-align: left; color: #aaa; }}
  td {{ padding: 5px 8px; border-bottom: 1px solid #223; }}
  tr:hover td {{ background: #1a2a4a; }}
  .log-box {{ background: #111; border-radius: 4px; padding: 12px; max-height: 240px; overflow-y: auto;
              font-size: 11px; line-height: 1.5; white-space: pre-wrap; word-break: break-all; }}
  .log-error {{ color: #f44336; }}
  .log-warn  {{ color: #ff9800; }}
  .log-info  {{ color: #90caf9; }}
  .ts {{ color: var(--gray); font-size: 11px; }}
  .full-width {{ grid-column: 1 / -1; }}
</style>
</head>
<body>
<header>
  <div>
    <h1>🧬 bioBakery Vis — Nextflow Monitor</h1>
    <div class="subtitle">Output: {output_dir} &nbsp;|&nbsp; Refreshes every 15 s &nbsp;|&nbsp; <span class="ts">{timestamp}</span></div>
  </div>
</header>

<div class="grid">

  <!-- Pipeline summary -->
  <div class="panel">
    <h2>Pipeline Status</h2>
    <div class="stat-row"><span>Total tasks</span><span class="badge badge-blue">{total}</span></div>
    <div class="stat-row"><span>Completed</span><span class="badge badge-green">{completed}</span></div>
    <div class="stat-row"><span>Running</span><span class="badge badge-yellow">{running}</span></div>
    <div class="stat-row"><span>Queued</span><span class="badge badge-gray">{queued}</span></div>
    <div class="stat-row"><span>Cached</span><span class="badge badge-blue">{cached}</span></div>
    <div class="stat-row"><span>Failed</span><span class="badge badge-red">{failed}</span></div>
    <div class="progress-bar">
      <div class="progress-fill" style="width:{progress_pct}%"></div>
    </div>
    <div class="ts" style="text-align:right">{progress_pct}% complete</div>
  </div>

  <!-- Host resources -->
  <div class="panel">
    <h2>Host Resources</h2>
    {resources_html}
  </div>

  <!-- SLURM jobs -->
  <div class="panel full-width">
    <h2>SLURM Queue (my jobs)</h2>
    {slurm_html}
  </div>

  <!-- Trace table -->
  <div class="panel full-width">
    <h2>Task Trace</h2>
    {trace_html}
  </div>

  <!-- Log tail -->
  <div class="panel full-width">
    <h2>Nextflow Log (last 60 lines)</h2>
    <div class="log-box">{log_html}</div>
  </div>

</div>
</body>
</html>
"""


def colorize_log(log_text: str) -> str:
    """Add simple color spans to log lines."""
    lines = []
    for line in log_text.splitlines():
        ll = line.lower()
        if "error" in ll or "failed" in ll or "exception" in ll:
            lines.append(f'<span class="log-error">{line}</span>')
        elif "warn" in ll:
            lines.append(f'<span class="log-warn">{line}</span>')
        elif "completed" in ll or "started" in ll or "finished" in ll:
            lines.append(f'<span class="log-info">{line}</span>')
        else:
            lines.append(line)
    return "\n".join(lines)


def render_slurm(jobs: list[dict]) -> str:
    if not jobs:
        return "<p style='color:#666'>SLURM not available or no jobs running.</p>"
    hdr = ["Job ID", "Name", "State", "Time", "Time Limit", "CPUs", "Memory", "Reason"]
    keys = ["job_id", "name", "state", "time", "time_limit", "cpus", "min_mem", "reason"]
    html = "<table><tr>" + "".join(f"<th>{h}</th>" for h in hdr) + "</tr>"
    for j in jobs:
        state = j.get("state", "")
        badge = ("badge-green" if state == "RUNNING"
                 else "badge-yellow" if state in ("PENDING", "COMPLETING")
                 else "badge-red")
        row = (f"<td>{j.get('job_id','')}</td>"
               f"<td>{j.get('name','')}</td>"
               f"<td><span class='badge {badge}'>{state}</span></td>"
               + "".join(f"<td>{j.get(k,'')}</td>" for k in keys[3:]))
        html += f"<tr>{row}</tr>"
    return html + "</table>"


def render_trace(rows: list[dict]) -> str:
    if not rows:
        return "<p style='color:#666'>No trace data yet. Run Nextflow with -with-trace.</p>"
    show_cols = ["task_id", "name", "status", "duration", "realtime",
                 "%cpu", "%mem", "rss", "peak_rss"]
    avail = [c for c in show_cols if any(c in r for r in rows)]
    html = "<table><tr>" + "".join(f"<th>{c}</th>" for c in avail) + "</tr>"
    for row in rows[-30:]:  # last 30 tasks
        status = row.get("status", "")
        badge = ("badge-green" if status == "COMPLETED"
                 else "badge-red" if status == "FAILED"
                 else "badge-yellow" if status in ("RUNNING", "SUBMITTED")
                 else "badge-blue" if status == "CACHED"
                 else "badge-gray")
        cells = ""
        for c in avail:
            val = row.get(c, "")
            if c == "status":
                cells += f"<td><span class='badge {badge}'>{val}</span></td>"
            else:
                cells += f"<td>{val}</td>"
        html += f"<tr>{cells}</tr>"
    return html + "</table>"


def render_resources(res: dict) -> str:
    if not HAS_PSUTIL:
        return "<p style='color:#666'>Install psutil for resource monitoring: pip install psutil</p>"
    cpu = res.get("cpu_pct", 0)
    mem = res.get("mem_pct", 0)
    html = (
        f'<div class="stat-row"><span>CPU usage</span>'
        f'<span><div class="progress-bar" style="width:120px;display:inline-block">'
        f'<div class="progress-fill" style="width:{cpu}%;background:{"#f44" if cpu>90 else "#ff9800" if cpu>70 else "#4caf50"}"></div>'
        f'</div> {cpu}%</span></div>'
        f'<div class="stat-row"><span>Memory usage</span>'
        f'<span><div class="progress-bar" style="width:120px;display:inline-block">'
        f'<div class="progress-fill" style="width:{mem}%;background:{"#f44" if mem>90 else "#ff9800" if mem>70 else "#4caf50"}"></div>'
        f'</div> {res.get("mem_used_gb","?")} / {res.get("mem_total_gb","?")} GB</span></div>'
    )
    return html


# ── FastAPI app ───────────────────────────────────────────────────────────────

def make_app(output_dir: Path, log_file: Path) -> "FastAPI":
    app = FastAPI(title="bioBakery Vis Monitor")

    trace_file = output_dir / "pipeline_info" / "trace.tsv"

    def _render_dashboard() -> str:
        trace_rows = read_trace(trace_file)
        summary = pipeline_summary(trace_rows)
        total = summary["total"] or 1
        progress = round((summary["completed"] + summary["cached"]) / total * 100)
        return DASHBOARD_HTML.format(
            output_dir=str(output_dir),
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total=summary["total"],
            completed=summary["completed"],
            running=summary["running"],
            queued=summary["queued"],
            cached=summary["cached"],
            failed=summary["failed"],
            progress_pct=progress,
            resources_html=render_resources(system_resources()),
            slurm_html=render_slurm(slurm_jobs()),
            trace_html=render_trace(trace_rows),
            log_html=colorize_log(read_log_tail(log_file)),
        )

    @app.get("/", response_class=HTMLResponse)
    def dashboard():
        return _render_dashboard()

    @app.get("/api/status")
    def api_status():
        trace_rows = read_trace(trace_file)
        return JSONResponse({
            "summary": pipeline_summary(trace_rows),
            "resources": system_resources(),
            "slurm_jobs": slurm_jobs(),
            "timestamp": datetime.now().isoformat(),
        })

    @app.get("/api/trace")
    def api_trace():
        return JSONResponse(read_trace(trace_file))

    @app.get("/api/log")
    def api_log():
        return {"log": read_log_tail(log_file)}

    return app


# ── CLI entry point ───────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="bioBakery Nextflow visualization workflow monitor"
    )
    parser.add_argument("--port", type=int, default=8050,
                        help="Dashboard port [default: 8050]")
    parser.add_argument("--output", default="vis_output",
                        help="Nextflow --output directory to watch [default: vis_output]")
    parser.add_argument("--log", default=".nextflow.log",
                        help="Nextflow log file [default: .nextflow.log]")
    parser.add_argument("--host", default="0.0.0.0",
                        help="Host to bind [default: 0.0.0.0]")
    args = parser.parse_args()

    if not HAS_FASTAPI:
        print("ERROR: FastAPI and uvicorn are required.")
        print("Install with: pip install fastapi uvicorn psutil")
        raise SystemExit(1)

    output_dir = Path(args.output).resolve()
    log_file = Path(args.log).resolve()

    print(f"bioBakery Vis Monitor")
    print(f"  Output dir : {output_dir}")
    print(f"  Log file   : {log_file}")
    print(f"  Dashboard  : http://localhost:{args.port}")
    print(f"  API        : http://localhost:{args.port}/api/status")
    print()

    app = make_app(output_dir, log_file)
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
