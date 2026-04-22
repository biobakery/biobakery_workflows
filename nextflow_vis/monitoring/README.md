# Monitoring Dashboard

Real-time web dashboard for bioBakery visualization Nextflow runs.

## Install

```bash
pip install fastapi uvicorn psutil
```

## Start

```bash
# From nextflow_vis/ directory:
python monitoring/dashboard.py \
    --output vis_output \
    --port   8050

# Open: http://localhost:8050
```

## Panels

| Panel | Source | Notes |
|-------|--------|-------|
| Pipeline status | `vis_output/pipeline_info/trace.tsv` | Tasks completed / running / failed |
| Host resources | `psutil` | CPU %, RAM GB — updates on refresh |
| SLURM queue | `squeue --me` | Only appears when running on SLURM |
| Task trace | `trace.tsv` | Last 30 tasks, sortable |
| Nextflow log | `.nextflow.log` | Last 60 lines, color-coded |

Auto-refreshes every 15 seconds (via `<meta http-equiv="refresh">`).

## JSON API

```bash
curl http://localhost:8050/api/status   # pipeline summary + SLURM + resources
curl http://localhost:8050/api/trace    # full trace table as JSON
curl http://localhost:8050/api/log      # last 60 log lines
```

## HPC port-forwarding

When the dashboard runs on an HPC login or head node:

```bash
# On your laptop:
ssh -L 8050:localhost:8050 user@hpc-cluster.edu

# Then open: http://localhost:8050
```
