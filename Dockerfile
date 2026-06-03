FROM python:3.12

ENV APP_HOME=/app
WORKDIR $APP_HOME

ENV PYTHONUNBUFFERED=1
ENV PORT=8080

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8080

# Listen on all interfaces; honor PORT from the host (e.g. Cloud Run, Knative).
# asgi:app is the Shiny app plus a /healthz readiness route (see asgi.py); it
# replaces `shiny run app.py` so the suite deploy contract's probe is served.
ENTRYPOINT ["sh", "-c", "exec uvicorn asgi:app --host 0.0.0.0 --port ${PORT:-8080}"]
