from django.core.urlresolvers import reverse
from django.shortcuts import render, get_object_or_404, redirect, get_list_or_404
from django.contrib import messages

# Create your views here.
from django.http import HttpResponseRedirect
from django.utils import timezone
from django.views.generic import CreateView

from .models import Dataset
from .forms import DatasetForm, EditDatasetForm


def index_view(request):
    # returns all measurements of the past, not of the future
    list_of_measurements = Dataset.objects.filter(measure_date__lte=timezone.now())#[:2]
    context = {
        'latest_measurement_list': list_of_measurements
    }
    return render(request, 'data_list.html', context)


def new_dataset(request):
    """
    Creates a new Dataset
    :param request:
    :return:
    """
    form = EditDatasetForm(request.POST or None, request.FILES or None)
    if form.is_valid():
        instance = form.save(commit=False)
        instance.save()
        messages.success(request, "Dataset saved successfully!")
        return HttpResponseRedirect(instance.get_absolute_url())
    context = {
        'form': form,
    }
    return render(request, "edit_dataset.html", context)


def edit_dataset(request, pk=None):
    """
    Updates a dataset.
    :param request:
    :param pk:
    :return:
    """
    instance = get_object_or_404(Dataset, pk=pk)
    form = EditDatasetForm(request.POST or None, request.FILES or None, instance=instance)
    if form.is_valid():
        instance = form.save(commit=False)
        instance.save()
        messages.success(request, "Dataset successfully updated!")
        return HttpResponseRedirect(instance.get_absolute_url())
    context = {
        'dataset': instance,
        'form': form,
    }
    return render(request, "edit_dataset.html", context)


def delete_dataset(request, pk=None):
    instance = get_object_or_404(Dataset, pk=pk)
    instance.delete()
    messages.success(request, "Dataset deleted successfully!")
    return redirect("index")


def detail_view(request, pk=None):
    """
    returns the detailed view of all the structures.
    """
    instance = get_object_or_404(Dataset, pk=pk)
    form = DatasetForm(instance=instance)
    context = {
        'dataset': instance,
        'form': form,
    }
    return render(request, 'detail.html', context)
















############################# old ###############
def list(request):
    # Handle file upload
    if request.method == 'POST':
        form = DatasetForm(request.POST, request.FILES)
        if form.is_valid():
            newdoc = Document(docfile=request.FILES['cif_file'])
            newdoc.save()
            # Redirect to the document list after POST
            return HttpResponseRedirect(reverse('stdb:list'))
    else:
        form = DatasetForm()  # A empty, unbound form

    # Load datasets and documents for the list page
    datasets = Dataset.objects.all()
    documents = Document.objects.all()
    #instance = get_object_or_404(Document, pk=?)
    context = {
        'datasets': datasets,
        'documents': documents,
        'form': form,
        #'instance': instance
    }
    # Render list page with the documents and the form
    return render(request, 'old_file_upload.html', context)
#########################################################