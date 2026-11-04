"""DeuringVUF SageMath implementation.

Protocol classes are loaded on demand so importing arithmetic helpers does not
initialize the public parameters or consume randomness.
"""

from importlib import import_module

__all__ = ["Keygen", "Eval", "Amortized_Eval", "Verify"]
_MODULES = {
    "Keygen": "protocol.keygen",
    "Eval": "protocol.eval",
    "Amortized_Eval": "protocol.eval",
    "Verify": "protocol.verify",
}


def __getattr__(name):
    if name not in _MODULES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{_MODULES[name]}", __name__), name)
    globals()[name] = value
    return value
