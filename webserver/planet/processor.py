from django.conf import settings

from natsort import natsorted #bette sorting for strings starting with numbers

import os
from os import path


def lst_experiments(request):

    request_full_path = request.get_full_path()

    lst_exp = []

    try:
        run_id = request_full_path.split("workspace/run/")[1].split("/")[0]
        exp_dir = path.join(settings.RUNS_DIR, run_id, "data", "experiments")
        lst_exp = list(natsorted(os.listdir(exp_dir)))
    except:
        pass


    context = {
      'lst_experiments': lst_exp
    }

    return context
