# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2016-05-14 21:10
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('stdb', '0004_auto_20160514_2057'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='refineresult',
            name='unit_cell',
        ),
    ]
