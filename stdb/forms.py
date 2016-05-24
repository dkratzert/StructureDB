from django import forms
#from .models import Dataset, Document


#from .models import Document
from stdb.models import Dataset



class DocumentForm(forms.ModelForm):
    cif_file = forms.FileField(label='Cif file', allow_empty_file=True, required=False)
    class Meta:
        model = Dataset
        fields = [str(i).split('.')[-1] for i in model._meta.fields]
        """
        fields = [
            'name',
            'flask_name',
            'machine',
            'measure_date',
            'formula',
            'operator',
            'cell_a',
            'cell_b',
            'cell_c',
            'alpha',
            'beta',
            'gamma',
            'cif_file',
        ]"""





"""
class CifDocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        #cif_file = forms.FileField(Document, label='Select a cif file')
        fields = [
            'cif_file',
        ]

class DocumentForm(forms.Form):
    docfile = forms.FileField(
        label='Select a file'
    )
"""