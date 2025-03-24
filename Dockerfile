FROM python:3.13

ENV APP_HOME /app
WORKDIR $APP_HOME

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True
ENV PORT 8080

COPY . .
RUN pip install --no-cache-dir --upgrade shiny
RUN pip install -r requirements.txt

# Run app on port 8080
EXPOSE ${PORT}
ENTRYPOINT ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "80"]

