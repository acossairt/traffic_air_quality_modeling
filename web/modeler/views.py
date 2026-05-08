import json
from pathlib import Path

from django.contrib import messages
from django.db import transaction
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_POST

from .forms import EDGE_SLIDER_FIELDS, EdgeFormSet, PatchFormSet, ProjectForm
from .models import Edge, Patch, Project, Run
from .tasks import start_run


def _network_payload(project):
    """Patches + edges in a JSON-friendly shape for the network preview."""
    return {
        "patches": [
            {"id": p.pk, "label": p.label, "initial_pop": p.initial_pop}
            for p in project.patches.all()
        ],
        "edges": [
            {"from": e.from_patch_id, "to": e.to_patch_id}
            for e in project.edges.all()
        ],
    }


def project_list(request):
    projects = Project.objects.all()
    return render(request, "modeler/project_list.html", {"projects": projects})


def project_create(request):
    if request.method == "POST":
        form = ProjectForm(request.POST)
        if form.is_valid():
            project = form.save()
            messages.success(request, f"Created project '{project.name}'.")
            return redirect("project_edit", pk=project.pk)
    else:
        form = ProjectForm()
    return render(request, "modeler/project_form.html", {"form": form, "project": None})


def project_edit(request, pk):
    """The main editor: project metadata + patch table + edge table."""

    project = get_object_or_404(Project, pk=pk)

    edge_kwargs = {"project": project}

    if request.method == "POST":
        if project.is_template:
            messages.error(request, "Templates are read-only — clone it first to make changes.")
            return redirect("project_edit", pk=project.pk)
        form = ProjectForm(request.POST, instance=project)
        patch_formset = PatchFormSet(request.POST, instance=project, prefix="patches")
        edge_formset = EdgeFormSet(
            request.POST, instance=project, prefix="edges", form_kwargs=edge_kwargs
        )

        if form.is_valid() and patch_formset.is_valid() and edge_formset.is_valid():
            form.save()
            patch_formset.save()
            edge_formset.save()
            messages.success(request, "Saved.")
            return redirect("project_edit", pk=project.pk)
    else:
        form = ProjectForm(instance=project)
        patch_formset = PatchFormSet(instance=project, prefix="patches")
        edge_formset = EdgeFormSet(
            instance=project, prefix="edges", form_kwargs=edge_kwargs
        )

    return render(
        request,
        "modeler/project_form.html",
        {
            "form": form,
            "project": project,
            "patch_formset": patch_formset,
            "edge_formset": edge_formset,
            "network_json": json.dumps(_network_payload(project)),
            "edge_slider_fields": EDGE_SLIDER_FIELDS,
        },
    )


def project_config_json(request, pk):
    """Inspectable: returns the config the Julia script will be fed."""
    project = get_object_or_404(Project, pk=pk)
    return JsonResponse(project.to_config())


@require_POST
def project_delete(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if project.is_template:
        messages.error(request, "Templates can't be deleted.")
        return redirect("project_list")
    name = project.name
    project.delete()
    messages.success(request, f"Deleted project '{name}'.")
    return redirect("project_list")


@require_POST
@transaction.atomic
def project_clone(request, pk):
    """Duplicate a project (typically the template) into a new editable sandbox."""
    source = get_object_or_404(Project, pk=pk)

    new_name = (request.POST.get("name") or "").strip() or f"Copy of {source.name}"

    clone = Project.objects.create(
        name=new_name,
        description=source.description,
        owner=request.user if request.user.is_authenticated else None,
        tspan_start=source.tspan_start,
        tspan_end=source.tspan_end,
        period=source.period,
        psi=source.psi,
        L=source.L,
        saveat=source.saveat,
        is_template=False,
    )

    patch_id_map = {}
    for p in source.patches.all():
        new_patch = Patch.objects.create(
            project=clone,
            label=p.label,
            initial_pop=p.initial_pop,
            display_order=p.display_order,
        )
        patch_id_map[p.pk] = new_patch

    for e in source.edges.all():
        Edge.objects.create(
            project=clone,
            from_patch=patch_id_map[e.from_patch_id],
            to_patch=patch_id_map[e.to_patch_id],
            vff=e.vff,
            nc_half_ff=e.nc_half_ff,
            vsharp=e.vsharp,
            demhf=e.demhf,
            shift=e.shift,
            dur=e.dur,
            dsharp=e.dsharp,
            bgd=e.bgd,
            onff=e.onff,
            on_half=e.on_half,
            onsharp=e.onsharp,
        )

    messages.success(request, f"Cloned '{source.name}' → '{clone.name}'.")
    return redirect("project_edit", pk=clone.pk)


@require_POST
def run_create(request, pk):
    project = get_object_or_404(Project, pk=pk)
    try:
        run = start_run(project)
    except ValueError as e:
        messages.error(request, str(e))
        return redirect("project_edit", pk=pk)
    return redirect("run_detail", pk=run.pk)


def run_detail(request, pk):
    run = get_object_or_404(Run, pk=pk)
    return render(request, "modeler/run_detail.html", {"run": run})


def run_status(request, pk):
    run = get_object_or_404(Run, pk=pk)
    return JsonResponse({
        "status": run.status,
        "error": run.error_message,
        "started_at": run.started_at.isoformat() if run.started_at else None,
        "completed_at": run.completed_at.isoformat() if run.completed_at else None,
    })


def run_results(request, pk):
    run = get_object_or_404(Run, pk=pk)
    if run.status != Run.STATUS_COMPLETED or not run.result_path:
        return JsonResponse({"error": "not ready", "status": run.status}, status=409)
    path = Path(run.result_path)
    if not path.exists():
        return JsonResponse({"error": "results file missing"}, status=410)
    with open(path) as f:
        return JsonResponse(json.load(f))
