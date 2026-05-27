from django import forms
from django.forms import inlineformset_factory

from .models import Edge, Patch, Project


# (min, max, step) ranges for each slider-controlled field on an Edge.
# Note: `shift` is presented to the user as a 24-hour-clock "peak time" with
# peak_time = shift + 6. The slider range below is in peak_time space.
EDGE_SLIDER_FIELDS = {
    "vff":        (0,    200,  1),
    "nc_half_ff": (0,    1,    0.01),
    "vsharp":     (0,    5,    0.1),
    "demhf":      (0,    20000, 100),
    "shift":      (0,    24,   0.1),
    "dur":        (-1,   1,    0.01),
    "dsharp":     (1,    500,  1),
    "bgd":        (0,    1,    0.01),
    "onff":       (0,    20000, 100),
    "on_half":    (0,    1,    0.01),
    "onsharp":    (0,    5,    0.1),
}

# Offset between the model's `shift` parameter and the user-facing "peak time".
# peak_time(24h) = shift + PEAK_TIME_OFFSET   (so shift=1 → 7 AM)
PEAK_TIME_OFFSET = 6


class ProjectForm(forms.ModelForm):
    class Meta:
        model = Project
        fields = [
            "name", "description",
            "tspan_start", "tspan_end",
            "period", "psi", "L", "saveat",
        ]
        widgets = {"description": forms.Textarea(attrs={"rows": 2})}


PatchFormSet = inlineformset_factory(
    Project,
    Patch,
    fields=["label", "initial_pop", "display_order"],
    extra=1,
    can_delete=True,
)


class EdgeForm(forms.ModelForm):
    class Meta:
        model = Edge
        fields = [
            "from_patch", "to_patch",
            "vff", "nc_half_ff", "vsharp",
            "demhf", "shift", "dur", "dsharp", "bgd",
            "onff", "on_half", "onsharp",
        ]

    def __init__(self, *args, project=None, **kwargs):
        super().__init__(*args, **kwargs)
        if project is not None:
            patches = project.patches.all()
            self.fields["from_patch"].queryset = patches
            self.fields["to_patch"].queryset = patches
        # Show only the patch label in the select; the default __str__ appends
        # "(project name)" which clutters the corridor card header.
        for name in ("from_patch", "to_patch"):
            if name in self.fields:
                self.fields[name].label_from_instance = lambda p: p.label

        # Tag the from-patch select so JS can identify it per card.
        if "from_patch" in self.fields:
            self.fields["from_patch"].widget.attrs["data-edge-from"] = "1"

        # Tag each numeric field so the live preview can read its value, and
        # set min/max/step so both the number input and its sibling slider use
        # consistent ranges.
        for name, (mn, mx, step) in EDGE_SLIDER_FIELDS.items():
            if name in self.fields:
                self.fields[name].widget.attrs.update({
                    "data-edge-param": name,
                    "min": mn,
                    "max": mx,
                    "step": step,
                })

        # User edits "peak time" (24h clock); we round-trip to the model's
        # `shift` parameter (peak_time = shift + PEAK_TIME_OFFSET). Override
        # data-edge-param after the loop above so the JS knows it's peak_time.
        if "shift" in self.fields:
            current_shift = (
                self.instance.shift
                if self.instance is not None and self.instance.pk is not None
                else (self.fields["shift"].initial or 0)
            )
            self.initial["shift"] = current_shift + PEAK_TIME_OFFSET
            self.fields["shift"].widget.attrs["data-edge-param"] = "peak_time"

    def clean_shift(self):
        # Convert back from peak_time to the underlying model `shift`.
        return self.cleaned_data["shift"] - PEAK_TIME_OFFSET


EdgeFormSet = inlineformset_factory(
    Project,
    Edge,
    form=EdgeForm,
    extra=1,
    can_delete=True,
)
