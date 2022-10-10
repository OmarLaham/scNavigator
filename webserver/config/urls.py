# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""

from django.conf import settings
from django.conf.urls import url
from django.urls import include, path, re_path
from django.conf.urls.static import static
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.contrib import admin
from django.views.generic import TemplateView
from django.views import defaults as default_views
#from rest_framework.authtoken.views import obtain_auth_token
from planet import views
from django.views.generic.base import RedirectView


urlpatterns = [
    path('admin/', admin.site.urls),          # Django admin route
    path('member/', include("apps.authentication.urls")), # Auth routes - login / register

    #main app URLs

    #homepage
    path("", views.workspace, name='home'),

    #start a run
    path("start", views.start, name="start"),

    #upload
    path("upload/<str:run_id>", views.upload, name="upload"),
    path("upload/<str:run_id>/push", views.workspace, name="upload_push"),

    #workspace
    path("workspace/run/<str:run_id>", views.workspace, name="workspace"),

    #saved degs lists & operations on them
    path("json_saved_degs_lists/<str:run_id>", views.workspace, name="saved_degs_lists"), #be careful: (s) at end in (name) param
    path("workspace/run/<str:run_id>/saved_degs_lists/<str:degs_list_id>", views.workspace, name="saved_degs_list"), #be careful: NO (s) at end in (name) param
    path("delete_degs_list/<str:run_id>/<str:degs_list_id>", views.workspace, name="delete_degs_list"),

    #running experiments
    path("run_epx_from_uploads/<str:run_id>/<str:upload_id>", views.workspace, name="run_exp_from_uploads"),
    path("run_exp_from_prev_exp/<str:run_id>/<str:prev_exp_id>/<str:cluster_id>", views.workspace, name="run_exp_from_prev_exp"),

    #r-scripts runners
    path("run_r_script/qc_metrics/<str:run_id>/<str:upload_name>/<int:min_nfeature_rna>/<int:max_nfeature_rna>/<int:percent_mt>", views.run_r_script_qc_metrics, name="r-qc-metrics"),
    path("run_r_script/elbow_plot/<str:run_id>/<str:upload_name>/<int:nDims>", views.workspace, name="r-elbow-plot"),

    #here you can pass a full saved degs list or any list of genes
    path("json_find_gene_intersection/<str:run_id>/<str:exp_id>/<str:cluster_id>/<str:genes_list>", views.workspace, name="find_gene_intersection")

    # UI Kits Html files
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) # server files from media

if settings.DEBUG:
    import debug_toolbar

    urlpatterns += [
        url(r'^__debug__/', include(debug_toolbar.urls)),
    ]

