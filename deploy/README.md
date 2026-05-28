# Deploy

Production deploy for `dimple-qc.odcambc.com` on Hetzner CPX21 (Ubuntu 24.04).

Topology: single Docker container behind Caddy on the same host. Caddy
terminates TLS, applies a 25 MB request body cap, applies a 60 req/min per-IP
rate limit, and reverse-proxies to the Shiny app on `127.0.0.1:8080`.

## First deploy

1. **Provision the server.** Hetzner Cloud Console → CPX21 → Ashburn → Ubuntu
   24.04 → SSH key. Note the IPv4 and IPv6.
2. **Add DNS records.** In whatever manages `odcambc.com`, add:
   ```
   dimple-qc.odcambc.com.  A     <server IPv4>
   dimple-qc.odcambc.com.  AAAA  <server IPv6>     # optional, recommended
   ```
   Verify with `dig +short dimple-qc.odcambc.com` from your laptop. Wait until
   it resolves before continuing — Caddy needs DNS to issue the TLS cert.
3. **SSH in and run setup.**
   ```bash
   ssh root@<server IPv4>
   curl -fL -o /tmp/setup.sh https://raw.githubusercontent.com/odcambc/dimple-qc-app/main/deploy/setup-server.sh
   bash /tmp/setup.sh
   ```
   The script is re-runnable. If it fails midway, re-run it.
4. **Verify.** From your laptop:
   ```bash
   curl -sSfI https://dimple-qc.odcambc.com | head -1   # expect HTTP/2 200
   ```

## Updating the app

```bash
ssh root@<server>
cd /opt/dimple-qc-app
git pull
docker build -t dimple-qc:latest .
docker rm -f dimple-qc
docker run -d --name dimple-qc --restart unless-stopped \
    -p 127.0.0.1:8080:8080 dimple-qc:latest
```

Caddy is not affected by app updates — it keeps its TLS cert and config.

## Operating

| Concern | Command |
|---|---|
| App logs | `docker logs --tail 100 -f dimple-qc` |
| Caddy logs | `journalctl -u caddy -n 100 -f` |
| Container status | `docker ps` |
| Restart app | `docker restart dimple-qc` |
| Restart Caddy | `systemctl restart caddy` |
| Reload Caddyfile after edit | `systemctl reload caddy` |
| Resource usage | `docker stats dimple-qc` |

## Session-isolation checklist (task #9 remaining items)

Run on the live server before announcing the URL:

1. **Single uvicorn worker.** Verify the container runs one Python process:
   ```bash
   docker exec dimple-qc ps -ef | grep python
   ```
   Should show one `shiny run app.py` process. If it shows multiple workers,
   session affinity is no longer guaranteed.
2. **Per-session tempdirs.** Upload a file in two browser sessions and verify:
   ```bash
   docker exec dimple-qc ls -la /tmp
   ```
   Each session should have its own directory under `/tmp/`. Confirm they are
   removed when the browser tab closes.
3. **Cache headers.** From your laptop:
   ```bash
   curl -sSI https://dimple-qc.odcambc.com | grep -i cache
   ```
   Expect `Cache-Control: private, no-store`.
4. **No file contents in logs.** Trigger a parse error (upload an empty
   file). Verify `docker logs dimple-qc` and `journalctl -u caddy` do not
   contain file contents or full upload paths.

## What's *not* set up here

- **Monitoring / alerting.** No Prometheus, no uptime check. For a small lab
  tool this is fine. If you want it, point a free UptimeRobot monitor at
  `https://dimple-qc.odcambc.com`.
- **Off-host backups.** Hetzner's weekly snapshot covers the host. Since the
  app is stateless (no database, no persistent uploads), there's nothing to
  back up beyond the host config.
- **Container auto-update.** Updates are manual. This is intentional — a
  user-facing tool shouldn't restart unexpectedly mid-session.
