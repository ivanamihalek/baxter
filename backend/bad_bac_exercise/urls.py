from django.urls import path, include
from rest_framework.routers import DefaultRouter
from bad_bac_exercise.views import AntibioticResMutationViewSet, DNASequenceView

router = DefaultRouter()
router.register(r'arm', AntibioticResMutationViewSet) # this works only for View sets

urlpatterns = [
    path('', include(router.urls)),
    path('fingerprint/', DNASequenceView.as_view()),  # Assuming it's a class-based view
]
