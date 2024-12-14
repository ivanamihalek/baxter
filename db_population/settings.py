import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# SECURITY WARNING: Modify this secret key if using in production!
SECRET_KEY = "6few3nci_q_o@l1dlbk81%wcxe!*6r29yu629&d97!hiqat9fa"

DEFAULT_AUTO_FIELD='django.db.models.AutoField'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'baxter',
        'USER': 'baxter',
        'PASSWORD': os.environ['MYSQL_PASSWORD'],
        'HOST': '127.0.0.1',   # an IP Address that your DB is hosted on
        'PORT': '3308',
    }
}

INSTALLED_APPS = ("models",)
