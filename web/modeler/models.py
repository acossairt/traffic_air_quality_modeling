from django.conf import settings
from django.db import models


class Project(models.Model):
    """A named scenario: a set of patches, edges, and simulation parameters."""

    name = models.CharField(max_length=200)
    description = models.TextField(blank=True)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        related_name="projects",
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    # simulation parameters (sim block in the JSON config)
    tspan_start = models.FloatField(default=0.0)
    tspan_end = models.FloatField(default=48.0)
    period = models.FloatField(default=24.0)
    psi = models.FloatField(default=100.0)
    L = models.FloatField(default=30.0)
    saveat = models.FloatField(default=0.1)

    class Meta:
        ordering = ["-updated_at"]

    def __str__(self):
        return self.name

    def to_config(self):
        """Build the JSON-serializable config that traffic_configurable.jl expects."""
        return {
            "patches": [p.to_config() for p in self.patches.all().order_by("display_order", "id")],
            "edges": [e.to_config() for e in self.edges.all().order_by("id")],
            "sim": {
                "tspan": [self.tspan_start, self.tspan_end],
                "period": self.period,
                "psi": self.psi,
                "L": self.L,
                "saveat": self.saveat,
            },
        }


class Patch(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name="patches")
    label = models.CharField(max_length=100)
    initial_pop = models.FloatField(default=0.0)
    display_order = models.IntegerField(default=0)

    class Meta:
        ordering = ["display_order", "id"]

    def __str__(self):
        return f"{self.label} ({self.project.name})"

    def to_config(self):
        # The Julia side uses integer ids; we expose Patch.pk directly so they
        # remain stable across edits.
        return {"id": self.pk, "label": self.label, "initial_pop": self.initial_pop}


class Edge(models.Model):
    """A directed corridor + its demand profile from one patch to another."""

    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name="edges")
    from_patch = models.ForeignKey(Patch, on_delete=models.CASCADE, related_name="outgoing")
    to_patch = models.ForeignKey(Patch, on_delete=models.CASCADE, related_name="incoming")

    # corridor params
    vff = models.FloatField(default=120.0, help_text="Free-flow velocity (km/h)")
    nc_half_ff = models.FloatField(default=0.3, help_text="Inverse capacity (corridor sat. fraction at half free-flow)")
    vsharp = models.FloatField(default=1.0, help_text="Sharpness of velocity-vs-density curve")

    # demand params
    demhf = models.FloatField(default=1000.0, help_text="Demand half-flux (saturation point in vehicles)")
    shift = models.FloatField(default=0.0, help_text="Phase shift of demand peak (hours)")
    dur = models.FloatField(default=0.97, help_text="Diurnal demand width")
    dsharp = models.FloatField(default=100.0, help_text="Sharpness of diurnal demand peak")
    bgd = models.FloatField(default=0.0, help_text="Background demand floor")

    # on-ramp / corridor entry params
    onff = models.FloatField(default=5000.0, help_text="On-ramp scale")
    on_half = models.FloatField(default=0.5, help_text="On-ramp half-saturation")
    onsharp = models.FloatField(default=1.0, help_text="On-ramp sharpness")

    class Meta:
        ordering = ["from_patch__display_order", "to_patch__display_order"]
        constraints = [
            models.UniqueConstraint(
                fields=["project", "from_patch", "to_patch"],
                name="unique_edge_per_pair",
            ),
        ]

    def __str__(self):
        return f"{self.from_patch.label} → {self.to_patch.label}"

    def to_config(self):
        return {
            "from": self.from_patch_id,
            "to": self.to_patch_id,
            "vff": self.vff,
            "nc_half_ff": self.nc_half_ff,
            "vsharp": self.vsharp,
            "demhf": self.demhf,
            "shift": self.shift,
            "dur": self.dur,
            "dsharp": self.dsharp,
            "bgd": self.bgd,
            "onff": self.onff,
            "on_half": self.on_half,
            "onsharp": self.onsharp,
        }


class Run(models.Model):
    STATUS_QUEUED = "queued"
    STATUS_RUNNING = "running"
    STATUS_COMPLETED = "completed"
    STATUS_FAILED = "failed"
    STATUS_CHOICES = [
        (STATUS_QUEUED, "Queued"),
        (STATUS_RUNNING, "Running"),
        (STATUS_COMPLETED, "Completed"),
        (STATUS_FAILED, "Failed"),
    ]

    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name="runs")
    status = models.CharField(max_length=16, choices=STATUS_CHOICES, default=STATUS_QUEUED)
    config_json = models.JSONField(help_text="Snapshot of the project config at run time")
    result_path = models.CharField(max_length=512, blank=True)
    error_message = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"Run {self.pk} ({self.project.name}, {self.status})"
