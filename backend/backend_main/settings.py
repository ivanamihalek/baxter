import os
# pip install django-configurations
# see https://django-configurations.readthedocs.io/en/stable/

from configurations import Configuration

class Base(Configuration):
    # all the base settings here...
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Quick-start development settings - unsuitable for production
    # SECURITY WARNING: keep the secret key used in production secret!
    SECRET_KEY = 'django-insecure-&exwl+fzxng80a&9f*&&b4r79jxfs0r@)!gai_=i!w@$gdgp(s'

    # Application definition
    INSTALLED_APPS = [
        'django.contrib.admin',
        'django.contrib.auth',
        'django.contrib.contenttypes',
        'django.contrib.sessions',
        'django.contrib.messages',
        'django.contrib.staticfiles',
        'django_extensions',
        'rest_framework',
        'corsheaders',
        'bad_bac_exercise'
    ]

    MIDDLEWARE = [
        'django.middleware.security.SecurityMiddleware',
        'django.contrib.sessions.middleware.SessionMiddleware',
        'django.middleware.common.CommonMiddleware',
        'django.middleware.csrf.CsrfViewMiddleware',
        'django.contrib.auth.middleware.AuthenticationMiddleware',
        'django.contrib.messages.middleware.MessageMiddleware',
        'django.middleware.clickjacking.XFrameOptionsMiddleware',
        'corsheaders.middleware.CorsMiddleware',
    ]


    ROOT_URLCONF = 'backend_main.urls'

    TEMPLATES = [
        {
            'BACKEND': 'django.template.backends.django.DjangoTemplates',
            'DIRS': [],
            'APP_DIRS': True,
            'OPTIONS': {
                'context_processors': [
                    'django.template.context_processors.debug',
                    'django.template.context_processors.request',
                    'django.contrib.auth.context_processors.auth',
                    'django.contrib.messages.context_processors.messages',
                ],
            },
        },
    ]

    WSGI_APPLICATION = 'backend_main.wsgi.application'

    # Password validation
    # https://docs.djangoproject.com/en/5.1/ref/settings/#auth-password-validators

    AUTH_PASSWORD_VALIDATORS = [
        {
            'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
        },
        {
            'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
        },
        {
            'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
        },
        {
            'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
        },
    ]

    # Internationalization
    # https://docs.djangoproject.com/en/5.1/topics/i18n/
    LANGUAGE_CODE = 'en-us'
    TIME_ZONE = 'UTC'
    USE_I18N = True
    USE_TZ = True

    # Static files (CSS, JavaScript, Images)
    # https://docs.djangoproject.com/en/5.1/howto/static-files/
    STATIC_URL = 'static/'

    # Default primary key field type
    # https://docs.djangoproject.com/en/5.1/ref/settings/#default-auto-field
    DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'


class Development(Base):
    # development settings here...
    DEBUG = True
    ALLOWED_HOSTS = ['localhost']
    CORS_ALLOWED_ORIGINS = [
        "http://localhost:3000",  # Next.js frontend
    ]
    # https://docs.djangoproject.com/en/5.1/ref/settings/#databases
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


class Production(Base):
    # production settings here...
    DEBUG = False
    ALLOWED_HOSTS = ['baxter-backend-ctnr']
    CORS_ALLOWED_ORIGINS = [
        "http://baxter-backend-ctnr:3000",  # Next.js frontend
    ]
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.mysql',
            'NAME': 'baxter',
            'USER': 'baxter',
            'PASSWORD': os.environ['MYSQL_PASSWORD'],
            'HOST': 'baxter-mariadb',  # Service name from docker-compose.prod.yml
            'PORT': '3306',  # Default MariaDB port

        }
    }
