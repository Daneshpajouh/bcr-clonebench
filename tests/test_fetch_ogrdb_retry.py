"""Tests for the download retry in bin/fetch_ogrdb.py, against a local stub HTTP server (no network)."""

import hashlib
import http.server
import io
import os
import sys
import threading
import unittest
from contextlib import redirect_stderr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import fetch_ogrdb as fo

BODY = b">IGHV1-2*02\nCAGGTG\n"


class Stub(http.server.BaseHTTPRequestHandler):
    plan = []  # status codes to send in turn; 200 sends BODY
    seen = 0

    def do_GET(self):
        code = self.plan[min(Stub.seen, len(self.plan) - 1)]
        Stub.seen += 1
        self.send_response(code)
        self.send_header("Content-Length", str(len(BODY) if code == 200 else 0))
        self.end_headers()
        if code == 200:
            self.wfile.write(BODY)

    def log_message(self, *args):
        pass


class Retry(unittest.TestCase):
    def setUp(self):
        Stub.seen = 0
        self.server = http.server.HTTPServer(("127.0.0.1", 0), Stub)
        threading.Thread(target=self.server.serve_forever, daemon=True).start()
        self.addCleanup(self.server.server_close)
        self.addCleanup(self.server.shutdown)
        self.url = f"http://127.0.0.1:{self.server.server_port}/set.fasta"
        self.waits = []

    def fetch(self, plan):
        Stub.plan = plan
        err = io.StringIO()
        with redirect_stderr(err):
            try:
                return fo.download(self.url, sleep=self.waits.append, timeout=10), err.getvalue()
            except ValueError as e:
                return e, err.getvalue()

    def test_two_500s_then_the_file(self):
        data, log = self.fetch([500, 500, 200])
        self.assertEqual(data, BODY)
        self.assertEqual((Stub.seen, self.waits), (3, [15, 30]))
        self.assertEqual(log.count("retrying in"), 2)
        self.assertIn("attempt 1 of 5 failed (HTTP 500", log)
        self.assertEqual(hashlib.md5(data).hexdigest(), hashlib.md5(BODY).hexdigest())

    def test_a_404_stops_at_once(self):
        err, log = self.fetch([404])
        self.assertIsInstance(err, ValueError)
        self.assertIn("HTTP 404", str(err))
        self.assertEqual((Stub.seen, self.waits, log), (1, [], ""))

    def test_five_500s_stop_after_the_last_wait(self):
        err, log = self.fetch([503])
        self.assertIsInstance(err, ValueError)
        self.assertIn("failed after 5 attempts; last: HTTP 503", str(err))
        self.assertEqual((Stub.seen, self.waits), (5, [15, 30, 60, 120]))

    def test_a_refused_connection_is_retried(self):
        self.server.shutdown()
        self.server.server_close()
        err = io.StringIO()
        with redirect_stderr(err):
            with self.assertRaises(ValueError) as ctx:
                fo.download(self.url, sleep=self.waits.append, timeout=10)
        self.assertEqual(self.waits, [15, 30, 60, 120])
        self.assertIn("URLError", str(ctx.exception))

    def test_the_digest_checks_still_stop_the_run(self):
        Stub.plan = [500, 200]
        row = {"md5": "md5:" + "0" * 32, "sha256": fo.NOT_PINNED}
        with redirect_stderr(io.StringIO()):
            data = fo.download(self.url, sleep=self.waits.append, timeout=10)
        with self.assertRaises(ValueError):
            fo.check_digests(data, row)


if __name__ == "__main__":
    unittest.main()
