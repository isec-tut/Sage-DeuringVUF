"""Parameter definitions, lazy initialization and random-state regressions."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import random
import subprocess
import sys
import json
import os

from src.config.parameters import DEFAULT_PARAMETERS


def test_scalar_parameters():
    assert DEFAULT_PARAMETERS.name == 'original'
    parameters = DEFAULT_PARAMETERS
    assert parameters.p == parameters.c*2**parameters.f-1
    assert parameters.N*parameters.odd_torsion_cofactor == parameters.p-1
    assert parameters.auxiliary_norm+parameters.N == 2**parameters.representation_exponent
    assert parameters.representation_cofactor*2**parameters.representation_exponent == parameters.torsion_order
    for change in ({'security_bits': 0}, {'c': 304}, {'f': 0},
                   {'odd_torsion_cofactor': 0}, {'odd_torsion_cofactor': 1183}):
        try:
            replace(parameters, **change)
        except ValueError:
            pass
        else:
            raise AssertionError(f'invalid parameters accepted: {change}')
    try:
        parameters.f = 239
    except FrozenInstanceError:
        pass
    else:
        raise AssertionError('parameter definitions must be immutable')


def test_lightweight_imports():
    script = '''
import random, sys
state = random.getstate()
from src.config import parameters
from src.config import constants
from src.config import context
from src.protocol import hashing
import precompute.parameters
assert random.getstate() == state
assert 'sage.all' not in sys.modules
assert 'precompute.build_context' not in sys.modules
assert constants.AUXILIARY_NORM == parameters.DEFAULT_PARAMETERS.auxiliary_norm
'''
    subprocess.run([sys.executable, '-c', script],
                   cwd=Path(__file__).resolve().parents[1], check=True)


def test_fixed_parameters():
    script = '''
import json, sys
from src.config.parameters import DEFAULT_PARAMETERS as p
from src.config import constants
from src.config import context
import precompute.parameters as legacy
assert (context.p, legacy.p) == (p.p, p.p)
assert constants.TORSION_EXPONENT == p.f
assert constants.REPRESENTATION_EXPONENT == p.representation_exponent
assert 'sage.all' not in sys.modules
print(json.dumps([p.name, p.c, p.f, p.secret_norm_strategy]))
'''
    # An obsolete selection variable must not change the fixed parameters.
    environment = dict(os.environ, DEURINGVUF_PARAMETERS='unsupported')
    result = subprocess.run([sys.executable, '-c', script],
                            cwd=Path(__file__).resolve().parents[1], env=environment,
                            capture_output=True, text=True, check=True)
    params = DEFAULT_PARAMETERS
    assert json.loads(result.stdout) == [params.name, params.c, params.f, params.secret_norm_strategy]


def test_context_initialization():
    from sage.all import set_random_seed
    from sage.misc.randstate import current_randstate
    from src.config import context
    from src.arithmetic.torsion import is_full_two_power_basis
    import precompute.parameters as legacy

    assert 'precompute.build_context' not in sys.modules
    random.seed(77)
    set_random_seed(77)
    python_state = random.getstate()
    sage_state = current_randstate()
    public = context.get_context()
    assert random.getstate() == python_state
    assert current_randstate() is sage_state
    assert context.get_context() is public
    from precompute.build_context import get_context as cached_context
    assert cached_context() is public
    assert cached_context(DEFAULT_PARAMETERS) is public
    assert legacy.O0 is context.O0 is public.O0
    assert context.U0_Fp4 is public.U0_Fp4
    assert public.parameters is DEFAULT_PARAMETERS
    assert public.alpha.reduced_norm() % public.N == 0
    assert public.I_P0.norm() == public.N
    assert is_full_two_power_basis(public.U0_Fp4, public.V0_Fp4, public.f)
    assert public.P0_Et.order() == public.Q0_Et.order() == public.N


def main():
    print(f'params: {DEFAULT_PARAMETERS.name}', flush=True)
    test_scalar_parameters()
    test_lightweight_imports()
    test_fixed_parameters()
    test_context_initialization()
    print('parameters: all tests passed')


if __name__ == '__main__':
    main()
