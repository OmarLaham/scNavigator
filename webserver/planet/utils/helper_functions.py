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

def str_to_numeric(string):
    value = None
    dtype = None

    try:
        value = int(string)
        dtype = "int"
        return value, dtype
    except:
        try:
            value = float(string)
            dtype = "float"
            return value, dtype
        except:
            return value, dtype

    return value, dtype

def is_valid_run(run_id):
    run_path = path.join(settings.RUNS_DIR, run_id)
    if path.exists(run_path):
        return True
    else:
        raise Exception("Internal Error: Invalid requested run id.")
