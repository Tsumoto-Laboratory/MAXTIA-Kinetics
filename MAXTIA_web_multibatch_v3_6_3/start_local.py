"""Serve the adjacent workbench on this computer and open its URL."""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import errno
import webbrowser


class LocalWorkbenchHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        if not self.path.startswith("/runtime/"):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()


def main():
    site = Path(__file__).resolve().parent / "dist"
    if not (site / "index.html").is_file():
        raise SystemExit("Extract the complete ZIP and keep start_local.py beside dist.")
    handler = partial(LocalWorkbenchHandler, directory=str(site))
    try:
        server = ThreadingHTTPServer(("127.0.0.1", 8002), handler)
    except OSError as exc:
        if exc.errno != errno.EADDRINUSE:
            raise
        server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    with server:
        url = f"http://127.0.0.1:{server.server_port}/?v=3.6.3"
        print(f"MAXTIA 3.6.3 · Open {url} in your browser. Press Ctrl+C to stop.", flush=True)
        webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
