from django.urls import path

from . import views

urlpatterns = [
    path("", views.project_list, name="project_list"),
    path("new/", views.project_create, name="project_create"),
    path("<int:pk>/", views.project_edit, name="project_edit"),
    path("<int:pk>/config.json", views.project_config_json, name="project_config_json"),
    path("<int:pk>/delete/", views.project_delete, name="project_delete"),
    path("<int:pk>/clone/", views.project_clone, name="project_clone"),
    path("<int:pk>/run/", views.run_create, name="run_create"),
    path("runs/<int:pk>/", views.run_detail, name="run_detail"),
    path("runs/<int:pk>/status", views.run_status, name="run_status"),
    path("runs/<int:pk>/results", views.run_results, name="run_results"),
]
