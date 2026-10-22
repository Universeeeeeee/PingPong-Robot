"""Download the pinned public sample and verify Git blob hashes."""
import concurrent.futures
import hashlib
import http.client
import json
from pathlib import Path
import urllib.request
import time

ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "source_manifest.json").read_text())

def fetch(entry):
    target = ROOT / "upstream" / entry["path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        existing = target.read_bytes()
        digest = hashlib.sha1(f"blob {len(existing)}\0".encode() + existing).hexdigest()
        if digest == entry["git_blob_sha1"]:
            return
    url = f"https://raw.githubusercontent.com/cogsys-tuebingen/spindoe/{manifest['commit']}/{entry['path']}"
    if entry["size"] > 2_000_000:
        def chunk(start):
            end = min(start + 1_000_000, entry["size"]) - 1
            part_path = ROOT / ".cache" / entry["git_blob_sha1"] / str(start)
            part_path.parent.mkdir(parents=True, exist_ok=True)
            if part_path.exists() and part_path.stat().st_size == end-start+1:
                return part_path.read_bytes()
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(request, timeout=120) as response:
                        if response.status != 206 or response.headers["Content-Range"] != f"bytes {start}-{end}/{entry['size']}":
                            raise ValueError("Server did not honor byte range")
                        part = response.read()
                    break
                except (OSError, http.client.HTTPException) as error:
                    if attempt == 3:
                        raise
                    print(f"Retry range {start}: {error}", flush=True)
                    time.sleep(1)
            if len(part) != end - start + 1:
                raise ValueError("Incomplete byte range")
            part_path.write_bytes(part)
            print(f"{entry['path']} bytes {start}-{end}", flush=True)
            return part
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            data = b"".join(pool.map(chunk, range(0, entry["size"], 1_000_000)))
    else:
        with urllib.request.urlopen(url, timeout=90) as response:
            data = response.read()
    digest = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    if len(data) != entry["size"] or digest != entry["git_blob_sha1"]:
        raise ValueError(f"Source verification failed: {entry['path']}")
    target.write_bytes(data)
    print(entry["path"], len(data), flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    list(pool.map(fetch, manifest["files"]))
