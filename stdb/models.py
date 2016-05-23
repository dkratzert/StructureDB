import datetime
import os

from django.core.urlresolvers import reverse
from django.core.validators import MinValueValidator, MaxValueValidator
from django.db import models

# Create your models here.
from django.db.models import Model
from django.forms import fields, forms
from django.http import HttpResponseRedirect
from django.shortcuts import render
from django.utils import timezone
from django.utils.html import escape
from django.utils.safestring import mark_safe


class Dataset(models.Model):
    name = models.CharField(max_length=200)
    flask_name = models.CharField(max_length=200)
    formula = models.CharField(max_length=200)
    machine = models.CharField(max_length=200)
    operator = models.CharField(max_length=200, verbose_name='Who measured the structure?')
    measure_date = models.DateTimeField(verbose_name='date measured')
    received = models.DateField(verbose_name='date received', blank=True, null=True)
    output = models.DateField(verbose_name='date outgoing', blank=True, null=True)
    crystal_size_x = models.FloatField(max_length=4, default=0, blank=True,
                                       validators = [MinValueValidator(0), MaxValueValidator(500)])
    crystal_size_y = models.FloatField(max_length=4, default=0, blank=True,
                                       validators=[MinValueValidator(0), MaxValueValidator(500)])
    crystal_size_z = models.FloatField(max_length=4, default=0, blank=True,
                                       validators=[MinValueValidator(0), MaxValueValidator(500)])
    customer = models.CharField(max_length=100, blank=True)
    is_publishable = models.BooleanField(default=False)
    comment = models.TextField(max_length=2000, blank=True)
    cell_a = models.FloatField(max_length=8, default=0, verbose_name='a', validators = [MinValueValidator(0), MaxValueValidator(500)])
    cell_b = models.FloatField(max_length=8, default=0, verbose_name='b', validators=[MinValueValidator(0), MaxValueValidator(500)])
    cell_c = models.FloatField(max_length=8, default=0, verbose_name='c', validators=[MinValueValidator(0), MaxValueValidator(500)])
    alpha = models.FloatField(max_length=8, default=0, validators = [MinValueValidator(0), MaxValueValidator(180)])
    beta = models.FloatField(max_length=8, default=0, validators = [MinValueValidator(0), MaxValueValidator(180)])
    gamma = models.FloatField(max_length=8, default=0, validators = [MinValueValidator(0), MaxValueValidator(180)])
    R1_all = models.FloatField(max_length=5, default=0,
                               validators = [MinValueValidator(0), MaxValueValidator(1)], blank=True)
    wR2_all = models.FloatField(max_length=5, default=0,
                                validators = [MinValueValidator(0), MaxValueValidator(1)], blank=True)
    R1_2s = models.FloatField(max_length=5, default=0,
                              validators = [MinValueValidator(0), MaxValueValidator(1)], blank=True)
    wR2_2s = models.FloatField(max_length=5, default=0,
                               validators = [MinValueValidator(0), MaxValueValidator(1)], blank=True)
    density = models.FloatField(max_length=5, verbose_name='density (calc)', blank=True, null=True)
    mu = models.FloatField(max_length=5, verbose_name='absorption [mm-1]', blank=True, null=True)
    formular_weight = models.FloatField(max_length=8, verbose_name='Formular weight', blank=True, null=True)
    colour = models.CharField(max_length=20, blank=True)
    shape = models.CharField(max_length=20, blank=True)
    temperature = models.FloatField(max_length=5, verbose_name='Temperature [K]', blank=True, default=0, null=True)
    crystal_system = models.CharField(max_length=15, blank=True)
    space_group = models.CharField(max_length=15, blank=True)
    volume = models.FloatField(max_length=5, blank=True, null=True)
    Z = models.IntegerField(blank=True, null=True)
    wavelength = models.FloatField(max_length=10, blank=True, null=True)
    radiation_type = models.CharField(max_length=18, blank=True)
    theta_min = models.FloatField(max_length=10, blank=True, null=True)
    theta_max = models.FloatField(max_length=10, blank=True, null=True)
    measured_refl = models.FloatField(max_length=10, blank=True, null=True)
    indep_refl = models.FloatField(max_length=20, blank=True, null=True)
    refl_used = models.FloatField(max_length=15, blank=True, null=True)
    r_int = models.FloatField(max_length=10, blank=True, null=True)
    parameters = models.FloatField(max_length=10, blank=True, null=True)
    restraints = models.FloatField(max_length=10, blank=True, null=True)
    peak = models.FloatField(max_length=10, blank=True, null=True)
    hole = models.FloatField(max_length=10, blank=True, null=True)
    goof = models.FloatField(max_length=10, blank=True, null=True)

    cell_a.short_description = 'Unit Cell Parameter a'
    cell_b.short_description = 'Unit Cell Parameter b'
    cell_c.short_description = 'Unit Cell Parameter c'
    alpha.short_description = 'Unit Cell Parameter alpha'
    beta.short_description = 'Unit Cell Parameter beta'
    gamma.short_description = 'Unit Cell Parameter gamma'
    Z.empty_value_display = '?'

    class Meta:
        ordering = ['-measure_date', '-name']

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('stdb:detail', kwargs={'pk': self.id})

    def was_measured_recently(self):
        now = timezone.now()
        return now - datetime.timedelta(days=7) <= self.measure_date <= now

    is_publishable.short_description = 'publishable?'
    is_publishable.boolean = True
    was_measured_recently.admin_order_field = 'measure_date'
    was_measured_recently.boolean = True
    was_measured_recently.short_description = 'Measured recently?'

    def publishable(self):
        if self.is_publishable:
            return True
        else:
            return False

    def format_formula(self):
        """
        Taken from: http: // www.mzan.com / article / 34067043 - django -
        not -rendering - html -
        with-tags - assigned - in -python - code.shtml  # sthash.QXxIeVLy.dpuf
        :param string: Chemical formula
        :return: formula as html string
        """
        new_string = ''
        for item in self.formula:
            if item.isdigit():
                # we know item is a number, so no need to escape it
                new_string += '<sub>'+ item +'</sub>'
            else:
                new_string += escape(item)
                # we built new_string from safe parts, so we can mark it as
                # safe to prevent autoescaping
        return mark_safe(new_string)


def get_filename(instance, name):
    """
    returns the path to where the uploaded file should be stored
    :param instance: An instance of the model where the FileField is defined.
    :param name: Filename to be stored
    :return: full path to store into
    """
    return os.path.join('cif_file/', os.path.join(instance.dataname.name, name))


class Document(models.Model):
    dataname = models.ForeignKey(Dataset, on_delete=models.CASCADE)#, null=True) # null=True alllows upload to no dataset!
    cif_file = models.FileField(upload_to=get_filename, verbose_name='CIF File', null=True, blank=True)
    res_file = models.FileField(upload_to=get_filename, verbose_name='RES File', null=True, blank=True)

    def __str__(self):
        return self.dataname.name

    def getsize(self):
        if 1000 <= self.cif_file.size < 1000000:
            return '{:.2f} {}'.format(self.cif_file.size / 1000, 'kB')
        else:
            return '{:.2f} {}'.format(self.cif_file.size / 1000000, 'MB')