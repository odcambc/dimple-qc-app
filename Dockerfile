FROM python:3.12

# uv: copy the standalone binaries from the official image — no pip bootstrap needed.
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

ENV APP_HOME=/app
WORKDIR $APP_HOME

ENV PYTHONUNBUFFERED=1
ENV PORT=8080
# Compile bytecode for faster cold starts; copy rather than hardlink since the
# uv cache and the target venv can live on different layers/filesystems.
ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy

# Install dependencies from the lockfile only. This layer is cached and rebuilt
# solely when pyproject.toml / uv.lock change, not on every source edit.
# --frozen: never re-resolve; fail if the lock is stale (CI catches drift).
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-dev --no-install-project

COPY . .

EXPOSE 8080

# asgi:app is the Shiny app plus a /healthz readiness route (see asgi.py); it
# replaces `shiny run app.py` so the suite deploy contract's probe is served.
# `uv run` executes inside the synced environment; --no-sync skips a redundant
# re-resolve at container start.
ENTRYPOINT ["sh", "-c", "exec uv run --no-sync uvicorn asgi:app --host 0.0.0.0 --port ${PORT:-8080}"]
