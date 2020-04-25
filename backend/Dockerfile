# Lightweight Python image.
FROM python:3.7-slim

# Copy local code to the container image.
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
RUN pip install -r requirements.txt

# Run the web service on container startup.
# 1 CPU core and 8 threads.
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 main:app
