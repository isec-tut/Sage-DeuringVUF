"""The DeuringVUF input map, shared by evaluation and verification."""

import hashlib
from src.config.parameters import N


def f_in(x):
    """Return the existing domain-separated SHA-256 scalar pair modulo N."""
    r = int.from_bytes(hashlib.sha256(b"0" + x).digest(), "big") % N
    s = int.from_bytes(hashlib.sha256(b"1" + x).digest(), "big") % N
    return r, s
