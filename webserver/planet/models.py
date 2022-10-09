# -*- encoding: utf-8 -*-

from django.db import models

class Run(models.Model):
    run_id = models.CharField(max_length=6, unique=True)
    description = models.CharField(max_length=100, null=False, blank=False)
    email = models.EmailField(max_length = 254)
    timestamp = models.DateTimeField(auto_now_add=True)

class Experiment(models.Model):
    exp_id = models.ForeignKey(Run, on_delete=models.CASCADE)
    description = models.CharField(max_length=100, null=False, blank=False)
    timestamp = models.DateTimeField(auto_now_add=True)
