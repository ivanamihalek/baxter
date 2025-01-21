import random

from django.db.models import OuterRef, Exists
from django.shortcuts import render

from rest_framework import viewsets
from rest_framework.response import Response
from rest_framework.views import APIView

from bad_bac_exercise import models
from bad_bac_exercise.models import AntibioticResMutation, PDBStructure, Decoy, Pdb2Mutation
from bad_bac_exercise.serializers import AntibioticResMutationSerializer, DNASequenceSerializer


class AntibioticResMutationViewSet(viewsets.ModelViewSet):
    queryset = AntibioticResMutation.objects.all()
    serializer_class = AntibioticResMutationSerializer


class DNASequenceView(APIView):

    def get_dna_sequences(self):
        """
        Fetches and returns a set of DNA sequences.
        """
        # Define the subquery to check for PDBStructure objects with dist_to_drug < 9
        pdb_subquery = Pdb2Mutation.objects.filter(
            antibio_res_mutation=OuterRef('pk'),
            dist_to_drug__lt=9
        )

        # Filter AntibioticResMutation objects based on the existence of such PDBStructure objects
        arm_entries = AntibioticResMutation.objects.filter(
            Exists(pdb_subquery),
            assemblies__isnull=False,
         ).distinct()

        if not arm_entries:
            raise Exception('No valid AntibioticResMutation entries found')

        # Select a random entry
        arm_entry = random.choice(list(arm_entries))

        # Ensure there are PDB entries for this mutation
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()

        if not pdb_entries:
            raise Exception('No PDB entries found for the selected mutation')

        # Fetch a random fingerprint from the assemblies
        assmb_entry = arm_entry.assemblies.first()
        if not assmb_entry:
            raise Exception('No assembly found for the selected mutation')

        fingerprints = assmb_entry.fingerprint_set.all()
        if not fingerprints:
            raise Exception('No fingerprints found for the selected assembly')

        fingerprint = random.choice(fingerprints).dna_seq

        # Fetch random decoys
        decoys = Decoy.objects.all()
        if not decoys:
            raise Exception('No decoys found')

        decoy1 = random.choice(decoys).dna_seq
        decoy2 = random.choice(decoys).dna_seq

        # Serialize the data
        data = {'fingerprint': fingerprint, 'decoy1': decoy1, 'decoy2': decoy2}
        serializer = DNASequenceSerializer(data=data)
        if serializer.is_valid():
            return serializer.data
        else:
            raise Exception('Serialization failed')

    def get(self, request):
        try:
            data = self.get_dna_sequences()
            return Response(data)
        except Exception as e:
            return Response({'error': str(e)}, status=500)
