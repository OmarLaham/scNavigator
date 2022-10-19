# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""
from django.conf import settings
from os import path
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, get_object_or_404, redirect
from django.template import loader
from django.http import HttpResponse, JsonResponse, Http404, HttpResponseRedirect
from django import template
from django.views import View
from planet.models import Run

from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt

import pandas as pd
import numpy as np
import math
import random

import re #regex
import json
import gzip
from glob import glob, escape

from natsort import natsorted #bette sorting for strings starting with numbers

import shutil

import os
from os import path

#to run bash scripts
import subprocess
from subprocess import Popen, PIPE

import planet.utils.helper_functions as hf

default_pval_thresh = 0.05

def index(request):
    return HttpResponseRedirect("/workspace/run/1")

def start(request):

    if request.method == 'POST':
        print("email:", request.POST["email"])
        print("description:", request.POST["description"])
        context = {}
        html_template = loader.get_template('home/start.html')
        return HttpResponse(html_template.render(context, request))
    else:
        context = {}
        html_template = loader.get_template('home/start.html')
        return HttpResponse(html_template.render(context, request))

def upload(request, run_id):

    context = {}
    html_template = loader.get_template('home/upload.html')
    return HttpResponse(html_template.render(context, request))

def lst_saved_deg_lists(request, run_id):


    lists_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists")
    if not path.exists(lists_path):
        raise Exception("Internal Error: Can't find", lists_path)

    saved_deg_lists = list(sorted(os.listdir(lists_path)))
    saved_deg_lists = [list.split(".csv")[0] for list in saved_deg_lists]

    context = {
        'saved_deg_lists': saved_deg_lists
    }


    print("> Done: lst_saved_deg_lists")

    return JsonResponse(context)


def save_deg_list(request, run_id, exp_title, cluster, deg_list_id):

    #copy from one to another
    source_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea",
                            "cluster_{0}_dea.csv".format(cluster))
    if not path.exists(source_path):
        raise Exception("Internal Error: Can't find source:", source_path)

    saved_list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists", "{0}.csv".format(deg_list_id))


    shutil.copy(source_path, saved_list_path)

    context = {}

    print("> Done: save_deg_list")

    return JsonResponse(context)


def del_saved_deg_list(request, run_id, deg_list_id):

    list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists", "{0}.csv".format(deg_list_id))
    if not path.exists(list_path):
        raise Exception("Internal Error: Can't find", list_path)

    os.remove(list_path)

    print("> del_saved_deg_list")

    return JsonResponse({})

def saved_deg_list(request, run_id, deg_list_id):

    saved_deg_list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists", "{0}.csv".format(deg_list_id))
    if not path.exists(saved_deg_list_path):
        raise Exception("Internal Error: Can't find", saved_deg_list_path)

    df_saved_deg_list = pd.read_csv(saved_deg_list_path)

    # get columns ready for intersection
    colnames = list(df_saved_deg_list.columns)
    colnames[0] = "gene"
    df_saved_deg_list.columns = colnames
    df_saved_deg_list = df_saved_deg_list[(df_saved_deg_list["p_val"] <= default_pval_thresh) | (pd.isna(df_saved_deg_list["p_val"]))] #NA condition is important for manually created lists

    gene_rows = df_saved_deg_list.values.tolist()

    tbl_html = "<tbody>"
    for index, row in df_saved_deg_list.iterrows():
        gene = row["gene"]
        avg_log2FC = row["avg_log2FC"]
        p_val = row["p_val"]
        p_val_adj = row["p_val_adj"]

        tbl_html += "<tr class='text-center'>"
        tbl_html += "<td><a target='_blank' href='https://www.genecards.org/cgi-bin/carddisp.pl?gene={0}' class='text-primary'>{1}</td>".format(gene, gene)
        tbl_html += "<td>{0}</td>".format(avg_log2FC)
        tbl_html += "<td>{0}</td>".format(p_val)
        tbl_html += "<td>{0}</td>".format(p_val_adj)
        tbl_html += "</tr>"
    tbl_html += "</tbody>"

    context = {
        "listID": deg_list_id,
        "tblHTML": tbl_html
    }

    html_template = loader.get_template('home/saved-deg-list.html')
    return HttpResponse(html_template.render(context, request))

    print("> saved_deg_list")

    return JsonResponse(context)

def create_deg_list_manually(request, run_id, deg_list_id, genes):

    saved_deg_list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists", "{0}.csv".format(deg_list_id))
    if path.exists(saved_deg_list_path):
        raise Exception("Internal Error: List ID already in Use.")

    genes = genes.split(",")

    header = ",p_val,avg_log2FC,pct.1,pct.2,p_val_adj\n"
    lines = []

    for gene in genes:
        lines.append("{0},NA,NA,NA,NA,NA".format(gene))

    lines = "\n".join(lines)

    file_content = header + lines

    with open(saved_deg_list_path, "w") as f:
        f.write(file_content)

    print("> Done: create_deg_list_manually")

    return JsonResponse({})


def json_find_cluster_to_list_gene_intersection(request, run_id, exp_title, cluster, saved_deg_list_id):
    dea_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")
    cluster_dea_path = path.join(dea_path, "cluster_{0}_dea.csv".format(cluster))

    df_cluster_dea = pd.read_csv(cluster_dea_path)
    # get columns ready for intersection
    colnames = list(df_cluster_dea.columns)
    colnames[0] = "gene"
    df_cluster_dea.columns = colnames
    df_cluster_dea = df_cluster_dea[df_cluster_dea["p_val"] <= default_pval_thresh]

    saved_deg_list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists",
                                    "{0}.csv".format(saved_deg_list_id))

    df_saved_deg_list = pd.read_csv(saved_deg_list_path)

    # get columns ready for intersection
    colnames = list(df_saved_deg_list.columns)
    colnames[0] = "gene"
    df_saved_deg_list.columns = colnames
    df_saved_deg_list = df_saved_deg_list[(df_saved_deg_list["p_val"] <= default_pval_thresh) | (pd.isna(df_saved_deg_list["p_val"]))] #NA is important for manually saved deg lists

    # let's do intersection ;)

    cluster_genes = df_cluster_dea["gene"].values.tolist()
    saved_deg_list_genes = df_saved_deg_list["gene"].values.tolist()

    gene_intersection = [gene for gene in cluster_genes if gene in saved_deg_list_genes]

    context = {
        "intersection": gene_intersection
    }

    print("> Done: json_find_cluster_to_list_gene_intersection")

    return JsonResponse(context)

def json_find_exp_to_list_gene_intersection(request, run_id, exp_title, saved_deg_list_id):

    dea_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")
    clusters_dea_paths = list(natsorted(os.listdir(dea_path)))

    # don't forget to remove clusters.csv ;)
    clusters = [cluster_dea_path.split("_dea.csv")[0].split("cluster_")[1] for cluster_dea_path in clusters_dea_paths if cluster_dea_path != "clusters.csv"]

    dict_intersection_result = dict()

    for cluster in clusters:

        df_cluster_dea = pd.read_csv(path.join(dea_path, "cluster_{0}_dea.csv".format(cluster)))
        # get columns ready for intersection
        colnames = list(df_cluster_dea.columns)
        colnames[0] = "gene"
        df_cluster_dea.columns = colnames
        df_cluster_dea = df_cluster_dea[df_cluster_dea["p_val"] <= default_pval_thresh]

        saved_deg_list_path = path.join(settings.RUNS_DIR, run_id, "saved_deg_lists",
                                        "{0}.csv".format(saved_deg_list_id))
        df_saved_deg_list = pd.read_csv(saved_deg_list_path)

        # get columns ready for intersection
        colnames = list(df_saved_deg_list.columns)
        colnames[0] = "gene"
        df_saved_deg_list.columns = colnames
        df_saved_deg_list = df_saved_deg_list[(df_saved_deg_list["p_val"] <= default_pval_thresh) | (pd.isna(df_saved_deg_list["p_val"]))] #NA is important for manually saved deg lists

        # let's do intersection ;)

        cluster_genes = df_cluster_dea["gene"].values.tolist()
        saved_deg_list_genes = df_saved_deg_list["gene"].values.tolist()

        gene_intersection = [gene for gene in cluster_genes if gene in saved_deg_list_genes]

        dict_intersection_result[cluster] = gene_intersection

    context = {
        "dict_intersection_result": dict_intersection_result
    }

    print("> Done: json_find_exp_to_list_gene_intersection")

    return JsonResponse(context)

def json_find_cluster_to_db_gene_intersection(request, run_id, exp_title, db_species, db_phase):

    dea_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")
    clusters_dea_paths = list(sorted(os.listdir(dea_path)))

    # don't forget to remove clusters.csv ;)
    clusters = [cluster_dea_path.split("_dea.csv")[0].split("cluster_")[1] for cluster_dea_path in clusters_dea_paths if
                cluster_dea_path != "clusters.csv"]

    df_db = pd.read_csv(path.join(settings.MEDIA_DIR, "marker_dbs", db_species, "{0}.csv".format(db_phase)))
    df_db = df_db[df_db["p_val"] <= default_pval_thresh] #this is very important to make comparison and intersection fare ;)
    # get columns ready for intersection
    colnames = list(df_db.columns)
    colnames[0] = "gene"
    df_db.columns = colnames

    dict_intersection_result = dict()

    for cluster in clusters:

        df_cluster_dea = pd.read_csv(path.join(dea_path, "cluster_{0}_dea.csv".format(cluster)))
        # get columns ready for intersection
        colnames = list(df_cluster_dea.columns)
        colnames[0] = "gene"
        df_cluster_dea.columns = colnames
        df_cluster_dea = df_cluster_dea[(df_cluster_dea["p_val"] <= default_pval_thresh) | pd.isna(df_cluster_dea["p_val"])] #NA is important for manually saved deg lists
        dea_cluster_genes = df_cluster_dea["gene"].values.tolist()

        dea_clusters = clusters # better code readability
        db_clusters = list(sorted(df_db["cluster"].value_counts().keys()))

        #store values to build df_intersections for this cluster
        intersections = []

        # my invention :P -> calc summation of pval and pval_adj to sort final result better (when we have many clusters with same n Genes) ;) You care only about sorting them then, not about any threshold.
        statistcs = []

        for db_cluster in db_clusters:
            df_db_cluster = df_db[df_db["cluster"] == db_cluster]
            db_cluster_genes = df_db_cluster["gene"].values.tolist()

            # sum p_val, p_val_adj (for better sorting later when we have many clusters with same nGenes)
            sum_p_val = df_db_cluster["p_val"].sum()
            sum_p_val_adj = df_db_cluster["p_val_adj"].sum()

            statistcs.append([db_cluster, sum_p_val, sum_p_val_adj])

            gene_intersection = list(set([gene for gene in dea_cluster_genes if gene in db_cluster_genes]))

            if len(gene_intersection) == 0:
                continue
            else:
                intersections.append([db_cluster, len(gene_intersection), ",".join(gene_intersection)])

        df_intersection = pd.DataFrame(intersections, columns=["cluster", "n_genes", "genes"])
        df_statistics = pd.DataFrame(statistcs, columns=["cluster", "sum_p_val", "sum_p_val_adj"])

        # group by cluster using max so we don't have multiple entries of same cluster
        df_intersection = df_intersection.groupby("cluster").max().reset_index()
        df_intersection.sort_values(by=["n_genes"], ascending=False, inplace=True)

        #e.g.: we have rows with nGenes: 3,3,3,2,1. This will return 3,2 and then we will return in the result every entry with 3 and 2 nGenes sorted by p_val_adj then p_val
        n_genes_to_select = list(sorted(df_intersection["n_genes"].values.tolist(), reverse=True))[:2]

        df_intersection = df_intersection[df_intersection["n_genes"].isin(n_genes_to_select)]
        df_intersection.reset_index(inplace=True)
        del df_intersection["index"]

        # add sum_p_val and sum_p_val_adj to sort results ;)
        df_intersection["sum_p_val"] = ""
        df_intersection["sum_p_val_adj"] = ""
        # fill with value from db_statistics
        for index, row in df_intersection.iterrows():
            inter_cluster = row["cluster"]
            df_intersection["sum_p_val"][index] = df_statistics[df_statistics["cluster"] == inter_cluster].iloc[0][
                "sum_p_val"]
            df_intersection["sum_p_val_adj"][index] = df_statistics[df_statistics["cluster"] == inter_cluster].iloc[0][
                "sum_p_val_adj"]

        df_intersection.sort_values(by=["n_genes", "sum_p_val_adj", "sum_p_val"], inplace=True,
                                    ascending=[False, True, True])

        cluster_intersection_result = df_intersection[["cluster", "n_genes", "genes"]].values.tolist()
        dict_intersection_result[cluster] = cluster_intersection_result

    context = {
        "dict_intersection_result": dict_intersection_result
    }

    print("> Done: json_find_cluster_to_db_gene_intersection")

    return JsonResponse(context)

def workspace(request, run_id):

    # validate run id
    hf.is_valid_run(run_id)

    raw_uploads_path = path.join(settings.RUNS_DIR, run_id, "data", "raw_uploads")

    raw_uploads = tuple(list(sorted((os.listdir(raw_uploads_path))))) # must convert to tuple (immutable) to be hashable -> return in context

    context = {
        "raw_uploads": raw_uploads
    }

    html_template = loader.get_template('home/workspace.html')

    print("> Done: workspace")

    return HttpResponse(html_template.render(context, request))


def run_r_script_qc_metrics(request, run_id, upload_name, min_nfeature_rna, max_nfeature_rna, percent_mt):

    timestamp = hf.generate_unix_timestamp() # create plots with timestamps to solve cashing problems

    subprocess.run(
        ["conda", "run", "-n", "single-cell", "Rscript",
         path.join(settings.SCRIPTS_DIR, "qc_metrics.R"),
         run_id, upload_name, str(min_nfeature_rna), str(max_nfeature_rna), str(percent_mt), str(timestamp)])



    context = {
        "img_src": "http://localhost:8000/media/" + path.join(settings.RUNS_DIR, run_id, "tmp", "qc_metrics_{0}.png".format(timestamp)).split("/media/")[1]
    }

    print("> Done: run_r_script_qc_metrics")

    return JsonResponse(context)

def run_r_script_elbow_plot(request, run_id, upload_name, min_nfeature_rna, max_nfeature_rna, percent_mt, n_dims):

    timestamp = hf.generate_unix_timestamp()  # create plots with timestamps to solve cashing problems

    subprocess.run(
        ["conda", "run", "-n", "single-cell", "Rscript",
         path.join(settings.SCRIPTS_DIR, "elbow_plot.R"),
         run_id, upload_name, str(min_nfeature_rna), str(max_nfeature_rna), str(percent_mt), str(n_dims), str(timestamp)])

    context = {
        "img_src": "http://localhost:8000/media/" + path.join(settings.RUNS_DIR, run_id, "tmp", "elbow_plot_{0}.png".format(timestamp)).split("/media/")[1]
    }

    print("> Done: run_r_script_elbow_plot")

    return JsonResponse(context)

def run_r_script_run_experiment(request, run_id, exp_title, upload_name, min_nfeature_rna, max_nfeature_rna, percent_mt, n_dims, clustering_res):

    timestamp = hf.generate_unix_timestamp()  # create plots with timestamps to solve cashing problems

    subprocess.run(
        ["conda", "run", "-n", "single-cell", "Rscript",
         path.join(settings.SCRIPTS_DIR, "run_experiment.R"),
         run_id, str(exp_title), upload_name, str(min_nfeature_rna), str(max_nfeature_rna), str(percent_mt), str(n_dims), str(clustering_res), str(timestamp)])

    # create df_exp_parameters to save exp parameters so we can load them any time the user loads an experiment later
    # BE CAREFUL: here we use data.rds (run_experiment result), not the original upload_name
    df_exp_parameters = pd.DataFrame([[run_id, str(exp_title), upload_name, str(min_nfeature_rna), str(max_nfeature_rna), str(percent_mt), str(n_dims), str(clustering_res)]],
             columns = ["run_id", "exp_title", "upload_name", "min_nfeature_rna", "max_nfeature_rna", "percent_mt", "n_dims", "clustering_res"])
    df_exp_parameters.to_csv(path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "exp_parameters.csv"))

    #clusters
    exp_dea_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")
    if not path.exists(exp_dea_path):
        raise Exception("Internal Error: Can't find {0}".format(exp_dea_path))
    clusters = os.listdir((exp_dea_path))
    clusters = [cluster for cluster in clusters if cluster.startswith("cluster_") and cluster.endswith("_dea.csv")] #you can use glob simply, too =)
    clusters = [cluster.split("_")[1].split("_")[0] for cluster in clusters]
    clusters = list(natsorted(clusters))

    context = {

        "umap_img_src": "http://localhost:8000/media/" +
                        path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "umap_{0}.png".format(timestamp)).split("/media/")[1],
        "tsne_img_src": "http://localhost:8000/media/" +
                        path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "tsne_{0}.png".format(timestamp)).split("/media/")[1],
        "trajectory_img_src": "http://localhost:8000/media/" +
                        path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "trajectory_pseudotime_plot_{0}.png".format(timestamp)).split("/media/")[1],
        # "umap_cell_embeddings": "http://localhost:8000/media/" + path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "umap_cell_embedings.png").split("/media/")[1],
        # "tsne_cell_embeddings": "http://localhost:8000/media/" + path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "tsne_cell_embeddings.png").split("/media/")[1],
        "data_rds_href": "http://localhost:8000/media/" +
                         path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "data.rds").split("/media/")[1],
        "clusters": clusters
    }

    print("> Done: run_r_script_run_experiment")

    return JsonResponse(context)

def run_r_script_list_clusters_dea_first(request, run_id, exp_title):

    # run Rscript to get clusters and DEGs for first cluster only if not run earlier for this exp.
    if not path.exists(path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea", "clusters.csv")) \
        or not path.exists(path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea", "cluster_0_dea.csv")):
        subprocess.run(
            ["conda", "run", "-n", "single-cell", "Rscript",
             path.join(settings.SCRIPTS_DIR, "list_clusters_dea_first.R"),
             run_id, str(exp_title)])

    #after Rscript you will have clusters.csv and cluster_i.csv(s) in /run_dir/run_id/data/experiments/exp_title/dea

    exp_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")

    df_clusters_path = path.join(exp_path, "clusters.csv")
    if not path.exists(df_clusters_path):
        raise Exception("Internal Error: can't find:", df_clusters_path)

    clusters = pd.read_csv(df_clusters_path)["levels(data)"].values.tolist()

    #generate tbody HTML for first cluster and store in cluster_1_degs_html
    df_cluster_degs_path = path.join(exp_path, "cluster_{0}_dea.csv".format(clusters[0]))
    if not path.exists(df_cluster_degs_path):
        raise Exception("Internal Error: can't find:", df_cluster_degs_path)

    df_cluster_degs = pd.read_csv(df_cluster_degs_path, index_col=False)
    colnames = list(df_cluster_degs.columns)
    colnames[0] = "gene"
    df_cluster_degs.columns = colnames
    df_cluster_degs = df_cluster_degs[["gene", "avg_log2FC", "p_val", "p_val_adj"]]
    df_cluster_degs = df_cluster_degs[df_cluster_degs["p_val"] <= default_pval_thresh]
    cluster_0_degs_as_list = df_cluster_degs.values.tolist()

    cluster_0_degs_as_str = []
    for i in range(len(cluster_0_degs_as_list)):
        cluster_0_degs_as_str.append(",".join([str(val) for val in cluster_0_degs_as_list[i]]) + "\n")
    cluster_0_degs_as_str = "".join(cluster_0_degs_as_str)

    context = {
        "clusters": clusters,
        "cluster_0_degs_as_list": cluster_0_degs_as_list,
        "cluster_0_degs_as_str": cluster_0_degs_as_str
    }

    print("> Done: run_r_script_list_clusters_dea_first")

    return JsonResponse(context)

def run_r_script_dea_cluster(request, run_id, exp_title, cluster):

    #run Rscript to get DEGs for this cluster only if not run earlier.
    if not path.exists(path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea", "cluster_{0}_dea.csv".format(cluster))):
        subprocess.run(
            ["conda", "run", "-n", "single-cell", "Rscript",
             path.join(settings.SCRIPTS_DIR, "dea_cluster.R"),
             run_id, str(exp_title), str(cluster)])

    #after Rscript you will have clusters.csv and cluster_i.csv(s) in /run_dir/run_id/data/experiments/exp_title/dea
    exp_dea_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "dea")

    #generate tbody HTML for first cluster and store in cluster_1_degs_html
    df_cluster_degs_path = path.join(exp_dea_path, "cluster_{0}_dea.csv".format(cluster))
    if not path.exists(df_cluster_degs_path):
        raise Exception("Internal Error: can't find:", df_cluster_degs_path)

    df_cluster_degs = pd.read_csv(df_cluster_degs_path, index_col=False)
    colnames = list(df_cluster_degs.columns)
    colnames[0] = "gene"
    df_cluster_degs.columns = colnames
    df_cluster_degs = df_cluster_degs[["gene", "avg_log2FC", "p_val", "p_val_adj"]]
    df_cluster_degs = df_cluster_degs[df_cluster_degs["p_val"] <= default_pval_thresh]


    cluster_degs_as_list = df_cluster_degs.values.tolist()

    cluster_degs_as_str = []
    for i in range(len(cluster_degs_as_list)):
        cluster_degs_as_str.append(",".join([str(val) for val in cluster_degs_as_list[i]]) + "\n")
    cluster_degs_as_str = "".join(cluster_degs_as_str)

    context = {
        "cluster_degs_as_list": cluster_degs_as_list,
        "cluster_degs_as_str": cluster_degs_as_str
    }

    print("> Done: run_r_script_dea_cluster")

    return JsonResponse(context)

def run_r_script_subset_cluster(request, run_id, upload_name, exp_title, cluster, min_nfeature_rna, max_nfeature_rna, percent_mt, n_dims, clustering_res):

    folder_name = "{0}_{1}_cluster{2}".format(upload_name, exp_title, cluster)

    subprocess.run(
        ["conda", "run", "-n", "single-cell", "Rscript",
         path.join(settings.SCRIPTS_DIR, "subset_cluster.R"),
         run_id, str(upload_name), str(exp_title), str(cluster), str(min_nfeature_rna), str(max_nfeature_rna), str(percent_mt), str(n_dims), str(clustering_res)])

    print("> Done: run_r_script_subset_cluster")

    return JsonResponse({})

def json_load_exp(request, run_id, exp_title):

    df_exp_params_path = path.join(settings.RUNS_DIR, run_id, "data", "experiments", exp_title, "exp_parameters.csv")
    if not path.exists(df_exp_params_path):
        raise Exception("Internal Error: Can't find", df_exp_params_path)

    df_exp_params = pd.read_csv(df_exp_params_path)
    row = df_exp_params.loc[0]

    dict_exp_params = {
        "run_id": str(row["run_id"]),
        "exp_title": str(row["exp_title"]),
        "upload_name": row["upload_name"],
        "min_nfeature_rna": str(row["min_nfeature_rna"]),
        "max_nfeature_rna": str(row["max_nfeature_rna"]),
        "percent_mt": str(row["percent_mt"]),
        "n_dims": str(row["n_dims"]),
        "clustering_res": str(row["clustering_res"])
    }

    context = {
        "dict_exp_parameters": dict_exp_params
    }

    print("> Done: json_load_exp")

    return JsonResponse(context)

def downloads(request):

    context = {}
    context['segment'] = 'downloads'

    html_template = loader.get_template('downloads.html')
    print("> Done: downloads")
    return HttpResponse(html_template.render(context, request))


#@login_required(login_url="/login/")
# def pages(request):
#     context = {}
#     # All resource paths end in .html.
#     # Pick out the html file name from the url. And load that template.
#     try:
#
#         load_template = request.path.split('/')[-1]
#         context['segment'] = load_template
#         html_template = loader.get_template( load_template )
#
#         return HttpResponse(html_template.render(context, request))
#
#     except template.TemplateDoesNotExist:
#
#         html_template = loader.get_template('not_used/templates/page-404.html')
#         return HttpResponse(html_template.render(context, request))
#
#     except:
#
#         html_template = loader.get_template('not_used/templates/page-500.html')
#         return HttpResponse(html_template.render(context, request))
