from django.conf.urls import url, include
from django.conf import settings

from stdb import admin
from .views import (
    index_view,
    detail_view,
    list,
    dataset_create,
    dataset_update,
    delete_dataset,
    upload)

app_name = 'stdb'
urlpatterns = [
    url(r'^$', index_view, name='index'),
    url(r'^(?P<pk>\d+)/$', detail_view, name='detail'),
    url(r'^(?P<pk>\d+)/edit/$', dataset_update, name='update'),
    url(r'^(?P<pk>\d+)/delete/$', delete_dataset, name='delete'),
    url(r'^(?P<pk>\d+)/upload/$', upload, name='upload'),
    url(r'^new/$', dataset_create, name='new'),
    url(r'^list/$', list, name='list'),
    ]

"""
if settings.DEBUG:
    urlpatterns += patterns('',
        url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {
            'document_root': settings.MEDIA_ROOT,
        }),
   )
"""