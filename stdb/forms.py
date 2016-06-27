import datetime
from django import forms
from django.forms import DateInput, CheckboxInput, Select, ModelChoiceField
from stdb.models import Dataset, Machines, Colours


class DateInput(forms.SelectDateWidget):
    """
    creates a date widget to select the date for
    e.g. the measurement
    """
    input_type = 'date'

    def __init__(self):
        super(DateInput, self).__init__()
        past_year = datetime.date.today().year - 15
        year_field = range(int(past_year), int(datetime.date.today().year) + 3)
        self.years = year_field


"""
class Meta:
        model = Dataset

def __init__(self, *args, **kwargs):
    super(IceCreamStoreUpdateForm, self).__init__(*args, **kwargs)
    self.fields['phone'].required = True
    self.fields['description'].required = True
"""


class DatasetForm(forms.ModelForm):
    class Meta:
        model = Dataset
        fields = [str(i).split('.')[-1] for i in model._meta.fields]
        is_publishable = forms.BooleanField(required=False)
        widgets = {
            'measure_date': DateInput,
            'received': DateInput,
            'output': DateInput,
            'is_publishable': CheckboxInput,
        }
        fields.extend(widgets.keys())



class EditDatasetForm(forms.ModelForm):
    class Meta:
        model = Dataset
        fields = ['name', 'flask_name', 'operator', 'machine', 'measure_date', 'received', 'output', 'picked_at',
                  'operator', 'colour', 'shape', 'is_publishable', 'service_structure',
                  'cif_file',
                  'fcf_file',
                  'res_file',
                  'raw_file',
                  'p4p_file',
                  'abs_file',
                  'eps_file',
                  'ls_file',
                  'pdf_file',
                  'checkcif_file',
                  'hkl_file',
                  'sfrm_file',
                  'cht_file',
                  'other_file1',
                  'other_file2',
                  'other_file3',
                  ]
        widgets = {
            'measure_date': DateInput,
            'received': DateInput,
            'output': DateInput,
            'is_publishable': CheckboxInput,
            #'machine': ModelChoiceField(queryset=Dataset.Machines_set.all()) #related_name=machines
        }


class MachinesForm(forms.Form):
    class Meta:
        model = Dataset
        machine = forms.ModelChoiceField(queryset=Machines.objects.all())