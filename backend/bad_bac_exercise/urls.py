from django.urls import path, include
from rest_framework.routers import DefaultRouter
from bad_bac_exercise.views import AntibioticResMutationViewSet

router = DefaultRouter()
router.register(r'arm', AntibioticResMutationViewSet)

urlpatterns = [
    path('', include(router.urls)),
]