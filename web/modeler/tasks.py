"""In-process threaded job runner: spawns Julia as a subprocess.

Trade-offs vs. a real queue (RQ/Celery):
- No extra services to run (no Redis, no worker process).
- Runs live in the Django process; if Django restarts mid-run, the Run row
  stays "running" until manually retried. Acceptable for a small lab tool.
"""
import json
import subprocess
import threading

from django.conf import settings
from django.utils import timezone

from .models import Run


def start_run(project):
    """Create a Run row for this project's current state and kick off Julia.

    Returns the Run instance immediately; the actual simulation happens in a
    background thread.
    """
    if project.patches.count() == 0:
        raise ValueError("Project has no patches.")

    config = project.to_config()
    run = Run.objects.create(
        project=project,
        status=Run.STATUS_QUEUED,
        config_json=config,
    )
    threading.Thread(target=_execute, args=(run.pk,), daemon=True).start()
    return run


def _execute(run_id):
    """Worker: write config to disk, invoke Julia, parse output, update Run."""
    from django.db import connection

    try:
        run = Run.objects.get(pk=run_id)
        run.status = Run.STATUS_RUNNING
        run.started_at = timezone.now()
        run.save(update_fields=["status", "started_at"])

        out_dir = settings.RUN_OUTPUT_DIR
        config_path = out_dir / f"run_{run_id}_config.json"
        result_path = out_dir / f"run_{run_id}_output.json"

        with open(config_path, "w") as f:
            json.dump(run.config_json, f)

        proc = subprocess.run(
            [str(settings.JULIA_BINARY), str(settings.JULIA_SCRIPT),
             str(config_path), str(result_path)],
            capture_output=True, text=True, timeout=600,
        )

        if proc.returncode != 0:
            run.status = Run.STATUS_FAILED
            run.error_message = (proc.stderr or proc.stdout or "Julia exited nonzero")[:5000]
        else:
            run.status = Run.STATUS_COMPLETED
            run.result_path = str(result_path)

        run.completed_at = timezone.now()
        run.save()

    except Exception as e:
        try:
            run = Run.objects.get(pk=run_id)
            run.status = Run.STATUS_FAILED
            run.error_message = repr(e)[:5000]
            run.completed_at = timezone.now()
            run.save()
        except Exception:
            pass
    finally:
        connection.close()
