#!/usr/bin/env bash
# One-time bootstrap for a fresh Hetzner CPX21 (Ubuntu 24.04 LTS) running dimple-qc-app.
# Run as root or via sudo. Re-runnable: each step is idempotent.

set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/odcambc/dimple-qc-app.git}"
INSTALL_DIR="${INSTALL_DIR:-/opt/dimple-qc-app}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() { printf '\n=== %s ===\n' "$*"; }

log "1/7 system update + base packages"
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get -y upgrade
apt-get -y install ca-certificates curl gnupg git ufw unattended-upgrades

log "2/7 firewall: SSH + HTTP + HTTPS only"
ufw default deny incoming
ufw default allow outgoing
ufw allow 22/tcp
ufw allow 80/tcp
ufw allow 443/tcp
ufw --force enable

log "3/7 unattended security upgrades"
dpkg-reconfigure -fnoninteractive unattended-upgrades

log "4/7 Docker"
if ! command -v docker >/dev/null; then
	curl -fsSL https://get.docker.com | sh
fi

log "5/7 Caddy with caddy-ratelimit plugin"
# Install the apt package first to get systemd unit, user, /etc/caddy, log dir.
if ! command -v caddy >/dev/null; then
	install -m 0755 -d /etc/apt/keyrings
	curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/gpg.key' \
		| gpg --dearmor -o /etc/apt/keyrings/caddy-stable.gpg
	curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/debian.deb.txt' \
		| tee /etc/apt/sources.list.d/caddy-stable.list >/dev/null
	apt-get update
	apt-get -y install caddy
fi

# Replace the stock binary with one built from caddyserver.com/download
# that includes caddy-ratelimit. The api/download endpoint streams a
# Linux/amd64 binary with the listed plugins.
log "5b/7 swapping Caddy binary for ratelimit-enabled build"
systemctl stop caddy || true
curl -fL -o /usr/bin/caddy.new \
	"https://caddyserver.com/api/download?os=linux&arch=amd64&p=github.com%2Fmholt%2Fcaddy-ratelimit"
chmod +x /usr/bin/caddy.new
mv /usr/bin/caddy.new /usr/bin/caddy
caddy version

log "6/7 install Caddyfile"
install -m 0644 "$SCRIPT_DIR/Caddyfile" /etc/caddy/Caddyfile
systemctl enable caddy
# Caddy will try to provision a Let's Encrypt cert on start. If DNS isn't
# resolving to this host yet, it will retry — HTTPS just won't work until DNS
# propagates.
systemctl restart caddy

log "7/7 clone + build + run the app container"
if [ ! -d "$INSTALL_DIR" ]; then
	git clone "$REPO_URL" "$INSTALL_DIR"
else
	git -C "$INSTALL_DIR" pull --ff-only
fi

cd "$INSTALL_DIR"
docker build -t dimple-qc:latest .
docker rm -f dimple-qc 2>/dev/null || true
# Bind to 127.0.0.1 only — Caddy is the public ingress. The container is
# unreachable from outside the host.
docker run -d \
	--name dimple-qc \
	--restart unless-stopped \
	-p 127.0.0.1:8080:8080 \
	dimple-qc:latest

cat <<EOF

Setup complete. Verify:
  1. DNS:  dig +short dimple-qc.odcambc.com   # should return this server's IP
  2. App:  curl -sSf http://127.0.0.1:8080 >/dev/null && echo "container OK"
  3. TLS:  curl -sSfI https://dimple-qc.odcambc.com | head -1
  4. Logs: journalctl -u caddy -n 50 --no-pager
           docker logs --tail 50 dimple-qc

To deploy an update later:
  cd $INSTALL_DIR && git pull && docker build -t dimple-qc:latest . \\
    && docker rm -f dimple-qc \\
    && docker run -d --name dimple-qc --restart unless-stopped \\
       -p 127.0.0.1:8080:8080 dimple-qc:latest
EOF
