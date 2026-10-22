"""Fetch the public Sony GCS supplementary MP4 with resumable byte ranges."""
import concurrent.futures
import hashlib
import http.client
from pathlib import Path
import time
import urllib.request

ROOT = Path(__file__).resolve().parent
URL = "https://media.githubusercontent.com/media/SonyResearch/ace_public/main/02_gaze_control_system.mp4"
SIZE = 93082306
SHA256 = "a5e43d1229feeea87373348e87a6315f707d9d02269bf8f182f9df9c79822a30"
CHUNK = 1_000_000
PARTS = ROOT / ".sony_parts"
PARTS.mkdir(exist_ok=True)
TARGET = ROOT / "sony_gcs_official.mp4"
START = time.monotonic()

def fetch(start):
    end = min(start + CHUNK, SIZE)-1
    part_path = PARTS / f"{start:09d}"
    if part_path.exists() and part_path.stat().st_size == end-start+1:
        return part_path
    request = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    for attempt in range(3):
        if time.monotonic()-START > 840:
            raise TimeoutError("14-minute download budget reached; completed parts preserved")
        try:
            with urllib.request.urlopen(request, timeout=65) as response:
                if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{SIZE}":
                    raise ValueError("Unexpected byte range response")
                data = response.read()
            if len(data) != end-start+1:
                raise http.client.IncompleteRead(data, end-start+1-len(data))
            part_path.write_bytes(data)
            print(f"downloaded {start}-{end}", flush=True)
            return part_path
        except (OSError, http.client.HTTPException) as error:
            print(f"retry {start} attempt {attempt+1}: {error}", flush=True)
            if attempt == 2:
                raise
            time.sleep(1)

if TARGET.exists() and TARGET.stat().st_size == SIZE and hashlib.sha256(TARGET.read_bytes()).hexdigest() == SHA256:
    print(f"Already verified: {TARGET}")
else:
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(fetch, range(0, SIZE, CHUNK)))
    digest = hashlib.sha256()
    staged = TARGET.with_suffix(".mp4.part")
    with staged.open("wb") as out:
        for start in range(0, SIZE, CHUNK):
            data = (PARTS / f"{start:09d}").read_bytes()
            digest.update(data)
            out.write(data)
    if staged.stat().st_size != SIZE or digest.hexdigest() != SHA256:
        raise ValueError("Full source verification failed")
    staged.replace(TARGET)
    print(f"VERIFIED {TARGET} SHA256={SHA256} elapsed={time.monotonic()-START:.1f}s", flush=True)
