#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys


def main():
    """Run administrative tasks."""
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'backend_main.settings')
    # not that 'Development' refers to the class defined in settings.py
    os.environ.setdefault('DJANGO_CONFIGURATION', 'Development')
    try:
        from configurations.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import configurations.managements. Check https://django-configurations.readthedocs.io"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    main()
