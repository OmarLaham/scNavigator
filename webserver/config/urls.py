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
    path("", views.homepage, name='home')

    # UI Kits Html files
]

if settings.DEBUG:
    import debug_toolbar

    urlpatterns += [
        url(r'^__debug__/', include(debug_toolbar.urls)),
    ]

