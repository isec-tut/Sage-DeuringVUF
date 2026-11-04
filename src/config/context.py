"""Lazy access to the default DeuringVUF arithmetic context.

Scalar values are available immediately. Curves, quaternions and torsion bases
are initialized together on first access and shared with precompute.parameters.
"""

from src.config.parameters import DEFAULT_PARAMETERS, c, f, N, p

_RUNTIME_NAMES = [
    "B", "i", "j", "k", "O0", "alpha", "I_P0", "kappa",
    "Fp2", "E0_fp2", "Et_fp2", "Fp4", "Fp2d", "zeta4", "E0_Fp4", "Et_Fp4",
    "iso_Et_Fp4_to_E0_Fp4", "iso_E0_Fp4_to_Et_Fp4",
    "R0", "R0_Fp4", "R0_E0", "alpha_coeffs", "kappa_coeffs",
    "P0_fp4", "Q0_fp4", "P0_Et_Fp4", "Q0_Et_Fp4", "P0_Et", "Q0_Et",
    "E0_fp2d", "U0_Fp4", "V0_Fp4", "Ms", "U0_Fp2", "V0_Fp2",
    "E0_Fp2_kummer", "UmV0_Fp2", "U0_Kummer", "V0_Kummer", "UmV0_Kummer",
]
__all__ = ["p", "f", "c", "N", "DEFAULT_PARAMETERS", "get_context", "main", *_RUNTIME_NAMES]


def get_context(parameters=DEFAULT_PARAMETERS):
    """Return the cached context for a frozen scalar parameter set."""
    from precompute.build_context import get_context as build
    return build(parameters)


def __getattr__(name):
    if name not in _RUNTIME_NAMES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(get_context(), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))


def main():
    """Inspect the default scalar parameters and initialized curve context."""
    from src.config.parameters import main as print_parameters
    print_parameters()
    ctx = get_context()
    for name in ("E0_fp2", "Et_fp2", "U0_Fp4", "V0_Fp4", "P0_Et", "Q0_Et"):
        value = getattr(ctx, name)
        print(f"{name}: {value}")
        if name in ("U0_Fp4", "V0_Fp4", "P0_Et", "Q0_Et"):
            print(f"order: {value.order()}")


if __name__ == "__main__":
    main()
