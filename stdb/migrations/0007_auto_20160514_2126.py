# -*- coding: utf-8 -*-
# Generated by Django 1.9.5 on 2016-05-14 21:26
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('stdb', '0006_refineresult_unit_cell'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='refineresult',
            name='unit_cell',
        ),
        migrations.AddField(
            model_name='refineresult',
            name='cell',
            field=models.OneToOneField(default=0, on_delete=django.db.models.deletion.CASCADE, to='stdb.UnitCell', verbose_name='Unit Cell'),
        ),
    ]
