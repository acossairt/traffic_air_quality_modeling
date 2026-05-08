"""Seed the database with the example 2-patch config for quick UI testing."""
import json
from pathlib import Path

from django.conf import settings
from django.core.management.base import BaseCommand

from modeler.models import Edge, Patch, Project


class Command(BaseCommand):
    help = "Import julia_model/example_config.json as a Project named 'Example 2-patch'."

    def handle(self, *args, **opts):
        config_path = Path(settings.BASE_DIR).parent / "julia_model" / "example_config.json"
        with open(config_path) as f:
            cfg = json.load(f)

        project, created = Project.objects.get_or_create(
            name="Example 2-patch",
            defaults=dict(
                description="Imported from julia_model/example_config.json",
                tspan_start=cfg["sim"]["tspan"][0],
                tspan_end=cfg["sim"]["tspan"][1],
                period=cfg["sim"]["period"],
                psi=cfg["sim"]["psi"],
                L=cfg["sim"]["L"],
                saveat=cfg["sim"]["saveat"],
            ),
        )
        if not created:
            project.patches.all().delete()
            project.edges.all().delete()

        # Map JSON ids to Patch instances so edges can resolve from/to.
        json_id_to_patch = {}
        for i, p in enumerate(cfg["patches"]):
            patch = Patch.objects.create(
                project=project,
                label=p["label"],
                initial_pop=p["initial_pop"],
                display_order=i,
            )
            json_id_to_patch[p["id"]] = patch

        for e in cfg["edges"]:
            Edge.objects.create(
                project=project,
                from_patch=json_id_to_patch[e["from"]],
                to_patch=json_id_to_patch[e["to"]],
                vff=e["vff"],
                nc_half_ff=e["nc_half_ff"],
                vsharp=e["vsharp"],
                demhf=e["demhf"],
                shift=e["shift"],
                dur=e["dur"],
                dsharp=e["dsharp"],
                bgd=e["bgd"],
                onff=e["onff"],
                on_half=e["on_half"],
                onsharp=e["onsharp"],
            )

        self.stdout.write(self.style.SUCCESS(
            f"{'Created' if created else 'Refreshed'} project '{project.name}' "
            f"with {project.patches.count()} patches and {project.edges.count()} edges."
        ))
