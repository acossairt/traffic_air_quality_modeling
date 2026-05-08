from django.contrib import admin

from .models import Edge, Patch, Project, Run


class PatchInline(admin.TabularInline):
    model = Patch
    extra = 0


class EdgeInline(admin.TabularInline):
    model = Edge
    extra = 0
    fk_name = "project"


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ("name", "owner", "updated_at")
    inlines = [PatchInline, EdgeInline]


@admin.register(Patch)
class PatchAdmin(admin.ModelAdmin):
    list_display = ("label", "project", "initial_pop", "display_order")
    list_filter = ("project",)


@admin.register(Edge)
class EdgeAdmin(admin.ModelAdmin):
    list_display = ("project", "from_patch", "to_patch", "vff", "demhf", "shift")
    list_filter = ("project",)


@admin.register(Run)
class RunAdmin(admin.ModelAdmin):
    list_display = ("id", "project", "status", "created_at", "completed_at")
    list_filter = ("status", "project")
    readonly_fields = ("config_json", "created_at", "started_at", "completed_at")
