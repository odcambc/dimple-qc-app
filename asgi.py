"""ASGI entrypoint: the Shiny Express app plus a /healthz readiness endpoint.

This exists to satisfy the suite-wide deploy contract (which requires every app
to expose an HTTP readiness probe). See docs/deploy-contract.md in the dms-tools
deploy repo.

`app.py` is a Shiny Express app, so it has no explicit `App` object to attach a
route to. `wrap_express_app()` builds that `App` for us. The returned `App` is
itself the ASGI entrypoint and owns the Shiny lifespan, so we keep it as the
root application and register one extra route on its internal Starlette router.
Mounting Shiny under a *separate* parent Starlette app would drop Shiny's
lifespan (Starlette does not propagate lifespan into mounted sub-apps).

Local dev is unchanged: `shiny run app.py` still works and simply omits
/healthz. The container launches this module instead: `uvicorn asgi:app`.
"""

from pathlib import Path

from shiny.express import wrap_express_app
from starlette.responses import JSONResponse
from starlette.routing import Route

app = wrap_express_app(Path(__file__).parent / "app.py")


async def healthz(_request):
    """Readiness probe. 200 once the app is serving; body is unspecified."""
    return JSONResponse({"status": "ok"})


# Insert ahead of Shiny's catch-all Mount("/") so /healthz resolves here and is
# not swallowed by Shiny's own routing.
app.starlette_app.router.routes.insert(
    0, Route("/healthz", healthz, methods=["GET"])
)
