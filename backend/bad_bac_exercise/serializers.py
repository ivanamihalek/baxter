import random

from rest_framework import serializers
from bad_bac_exercise.models import AntibioticResMutation, PDBStructure, Decoy


class AntibioticResMutationSerializer(serializers.ModelSerializer):
    class Meta:
        model = AntibioticResMutation
        fields = ['id', 'mutation', 'gene']


class DNASequenceSerializer(serializers.Serializer):
    """Serializer for DNA sequences."""
    fingerprint = serializers.CharField()
    decoy1 = serializers.CharField()
    decoy2 = serializers.CharField()

    class Meta:
        fields = ['fingerprint', 'decoy1', 'decoy2']
