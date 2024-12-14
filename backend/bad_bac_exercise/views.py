from django.shortcuts import render

from rest_framework import viewsets
from bad_bac_exercise.models import AntibioticResMutation
from bad_bac_exercise.serializers import AntibioticResMutationSerializer


class AntibioticResMutationViewSet(viewsets.ModelViewSet):
    queryset = AntibioticResMutation.objects.all()
    serializer_class = AntibioticResMutationSerializer
