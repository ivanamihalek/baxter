# Use the official Python image from the Docker Hub
FROM python:3.12-slim

# Set environment variables to prevent Python from writing .pyc files and to buffer stdout/stderr
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Create a directory for the application code
RUN mkdir /app

# Create a directory for the logging output
ENV LOGGING_DIR /scratch/logging
RUN mkdir -p /scratch/logging

# Set the working directory inside the container
WORKDIR /app

# Copy the requirements file into the container
COPY requirements.txt /app/

# will need those for the mysqlcline
RUN apt-get update -y
RUN apt-get install -y python3-dev build-essential
RUN apt-get install -y default-libmysqlclient-dev pkg-config
# Install the dependencies specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


# Copy the entire project directory into the container
COPY . /app/

# Expose port 8000 for the Django application
EXPOSE 8000

# this should come from CLI:  docker build --build-arg mysql_password=ajfghkaf Dockerfile
ARG mysql_password="nopass"
ENV MYSQL_PASSWORD=${mysql_password}
# Command to run the Django development server
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]

