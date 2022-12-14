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
    path("", views.index, name='home'),

    #start a run
    path("start", views.start, name="start"),

    #upload
    path("upload/<str:run_id>", views.upload, name="upload"),
    path("upload/<str:run_id>/push", views.workspace, name="upload_push"),

    #workspace
    path("workspace/run/<str:run_id>", views.workspace, name="workspace"),

    #saved degs lists & operations on them
    path("json_lst_saved_deg_lists/<str:run_id>", views.lst_saved_deg_lists, name="lst_saved_deg_lists"), #be careful: (s) at end in (name) param
    path("save_deg_list/<str:run_id>/<str:exp_title>/<str:cluster>/<str:deg_list_id>", views.save_deg_list, name="save_deg_list"), #be careful: NO (s) at end in (name) param
    path("del_saved_deg_list/<str:run_id>/<str:deg_list_id>", views.del_saved_deg_list, name="del_saved_deg_list"),
    path("saved_deg_list/<str:run_id>/<str:deg_list_id>", views.saved_deg_list, name="saved_deg_list"),
    path("create_deg_list_manually/<str:run_id>/<str:deg_list_id>/<str:genes>", views.create_deg_list_manually, name="create_deg_list_manually"),
    # #this view can be used for specific or all clusters if <str:clusters> == "all"
    # path("cross_with_saved_deg_lists/<str:run_id>/<str:exp_title>/<str:clusters>/<str:saved_deg_lists>", views.cross_with_saved_deg_lists, name="cross_with_saved_deg_lists"),


    #running experiments
    path("run_epx_from_uploads/<str:run_id>/<str:upload_id>", views.workspace, name="run_exp_from_uploads"),
    path("run_exp_from_prev_exp/<str:run_id>/<str:prev_exp_id>/<str:cluster_id>", views.workspace, name="run_exp_from_prev_exp"),

    #r-scripts runners
    path("run_r_script/qc_metrics/<str:run_id>/<str:upload_name>/<int:min_nfeature_rna>/<int:max_nfeature_rna>/<int:percent_mt>", views.run_r_script_qc_metrics, name="r-qc-metrics"),
    path("run_r_script/elbow_plot/<str:run_id>/<str:upload_name>/<int:min_nfeature_rna>/<int:max_nfeature_rna>/<int:percent_mt>/<int:n_dims>", views.run_r_script_elbow_plot, name="r-elbow-plot"),
    path("run_r_script/run_experiment/<str:run_id>/<str:exp_title>/<str:upload_name>/<int:min_nfeature_rna>/<int:max_nfeature_rna>/<int:percent_mt>/<int:n_dims>/<str:clustering_res>", views.run_r_script_run_experiment, name="r-run-exp"),
    path("run_r_script/list_clusters_dea_first/<str:run_id>/<str:exp_title>", views.run_r_script_list_clusters_dea_first, name="r-list-clusters-dea-first"),
    path("run_r_script/dea_cluster/<str:run_id>/<str:exp_title>/<str:cluster>", views.run_r_script_dea_cluster, name="r-run-exp"),
    #start experiment from cluster (subsetting)
    path("run_r_script/subset_cluster/<str:run_id>/<str:upload_name>/<str:exp_title>/<str:cluster>/<int:min_nfeature_rna>/<int:max_nfeature_rna>/<int:percent_mt>/<int:n_dims>/<str:clustering_res>", views.run_r_script_subset_cluster, name="r-subset-cluster"),

    #here you can pass a full saved degs list or any list of genes
    path("json_find_cluster_to_list_gene_intersection/<str:run_id>/<str:exp_title>/<str:cluster>/<str:saved_deg_list_id>", views.json_find_cluster_to_list_gene_intersection, name="json_find_cluster_to_list_gene_intersection"),
    path("json_find_exp_to_list_gene_intersection/<str:run_id>/<str:exp_title>/<str:saved_deg_list_id>", views.json_find_exp_to_list_gene_intersection, name="json_find_exp_to_list_gene_intersection"),
    path("json_find_cluster_to_db_gene_intersection/<str:run_id>/<str:exp_title>/<str:db_species>/<str:db_phase>", views.json_find_cluster_to_db_gene_intersection, name="json_find_cluster_to_db_gene_intersection"),

    #fetch current run experiments names as a list
    path("json_load_exp/<str:run_id>/<str:exp_title>", views.json_load_exp, name="json_load_exp")

    # UI Kits Html files
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) # server files from media

if settings.DEBUG:
    import debug_toolbar

    urlpatterns += [
        url(r'^__debug__/', include(debug_toolbar.urls)),
    ]

