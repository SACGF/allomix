"""Process-stable seed derivation for the paper's simulation scripts.

Python's built-in hash() is randomised per process for str/bytes (PEP 456), so
``seed + hash(str(frac))`` gave different simulated samples on every run. This
helper hashes the repr with SHA-256 instead, so a given (seed, fraction) always
maps to the same sample seed.
"""

import hashlib


def derive_seed(*parts: object) -> int:
    """Deterministic 32-bit seed from arbitrary parts (stable across processes)."""
    digest = hashlib.sha256(repr(parts).encode("utf-8")).digest()[:4]
    return int.from_bytes(digest, "big")
