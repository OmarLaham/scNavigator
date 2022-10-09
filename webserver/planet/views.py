# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""
from django.conf import settings
from os import path
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, get_object_or_404, redirect
from django.template import loader
from django.http import HttpResponse, JsonResponse, Http404
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

import os
from os import path

import planet.utils.helper_functions as hf



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

def workspace(request, run_id):

    # validate run id
    hf.is_valid_run(run_id)

    raw_uploads_path = path.join(settings.RUNS_DIR, run_id, "data", "raw_uploads")

    raw_uploads = tuple(list(sorted((os.listdir(raw_uploads_path))))) # must convert to tuple (immutable) to be hashable -> return in context

    context = {
        "raw_uploads": raw_uploads
    }

    html_template = loader.get_template('home/workspace.html')
    return HttpResponse(html_template.render(context, request))





def downloads(request):

    context = {}
    context['segment'] = 'downloads'

    html_template = loader.get_template('downloads.html')
    return HttpResponse(html_template.render(context, request))

def HomeJSONView(request):

    context = {}
    return JsonResponse(context)


class MetadataView(View):
    template_name = 'metadata.html'

    def get(self, request, *args, **kwargs):

        context = {}


        #context["full_data_table_thead_content"] = thead_content
        #context["full_data_table_tbody_content"] = tbody_content

        return render(request, self.template_name, context)

#@login_required(login_url="/login/")
def pages(request):
    context = {}
    # All resource paths end in .html.
    # Pick out the html file name from the url. And load that template.
    try:

        load_template = request.path.split('/')[-1]
        context['segment'] = load_template
        html_template = loader.get_template( load_template )

        return HttpResponse(html_template.render(context, request))

    except template.TemplateDoesNotExist:

        html_template = loader.get_template('not_used/templates/page-404.html')
        return HttpResponse(html_template.render(context, request))

    except:

        html_template = loader.get_template('not_used/templates/page-500.html')
        return HttpResponse(html_template.render(context, request))
