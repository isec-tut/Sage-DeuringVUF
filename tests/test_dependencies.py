"""Dependency integrity and connections to the local SQIsign fork.

Run: sage -python -m tests.test_dependencies
"""

import hashlib
import json
from pathlib import Path

from src.config import context
from src.arithmetic.kummer import KummerPoint
from src.protocol.keygen import Keygen
from src.protocol.eval import Eval, Amortized_Eval
from src.protocol.verify import Verify


def test_upstream_integrity():
    root = Path(__file__).resolve().parents[1] / 'third_party/sqisign_prism_v2'
    manifest = json.loads((root / 'UPSTREAM.json').read_text())
    modifications = manifest.get('local_modifications', {})
    assert set(modifications) <= set(manifest['files'])
    for name, expected in manifest['files'].items():
        if name in modifications:
            expected = modifications[name]['git_blob']
        content = (root / name).read_bytes()
        blob = b'blob ' + str(len(content)).encode() + b'\0' + content
        assert hashlib.sha1(blob).hexdigest() == expected, name



def test_module_boundaries():
    import src
    from src.arithmetic import torsion
    from third_party.sqisign_prism_v2 import hd
    from third_party.sqisign_prism_v2.theta.montgomery_isogenies.kummer_line import KummerPoint as ForkKummerPoint
    from third_party.sqisign_prism_v2.theta.theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
    from precompute import torsion_points

    assert context.U0_Kummer.__class__ is KummerPoint
    assert KummerPoint is ForkKummerPoint
    assert hd.EllipticProductIsogenySqrt is EllipticProductIsogenySqrt
    assert torsion.TorsionBasis_2e is torsion_points.TorsionBasis_2e
    assert (src.Keygen, src.Eval, src.Amortized_Eval, src.Verify) == (
        Keygen, Eval, Amortized_Eval, Verify,
    )


def main():
    test_upstream_integrity()
    test_module_boundaries()
    print('dependencies: all tests passed', flush=True)


if __name__ == '__main__':
    main()
