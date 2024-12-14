from rest_framework import serializers
from bad_bac_exercise.models import AntibioticResMutation


class AntibioticResMutationSerializer(serializers.ModelSerializer):
    class Meta:
        model = AntibioticResMutation
        fields = ['id', 'mutation', 'gene']
