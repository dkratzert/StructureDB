from django.contrib import admin

# Register your models here.
from .models import Dataset


class DatasetAdmin(admin.ModelAdmin):
    fieldsets = [
        ('Measurement', {'fields': ['name', 'measure_date', 'operator', 'flask_name', 'machine', 'received', 'output']}),
        ('Misc', {'fields': ['is_publishable', 'formula', 'Z', 'comment']}),
        ('Results', {'fields': [ 'cell_a', 'cell_b', 'cell_c',
                                'alpha', 'beta', 'gamma', 'R1_all', 'wR2_all', 'R1_2s', 'wR2_2s',
                                'density', 'mu', 'formular_weight', 'colour', 'shape', 'temperature', 'crystal_system',
                                'space_group', 'volume', 'wavelength', 'radiation_type', 'theta_min', 'theta_max',
                                'measured_refl', 'indep_refl', 'refl_used', 'r_int', 'parameters', 'restraints',
                                'peak', 'hole', 'goof'
                                ]
                     }
         )
    ]
    list_display = ('name', 'measure_date', 'was_measured_recently', 'is_publishable', 'machine', 'operator')
    list_filter = ['measure_date']
    search_fields = ['name', 'operator', 'flask_name']

admin.site.register(Dataset, DatasetAdmin)


"""
Add possibility to add machines like Choices in https://docs.djangoproject.com/en/1.9/intro/tutorial07/
"""
