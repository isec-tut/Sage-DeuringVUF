"""Run each test module in a fresh process: sage -python -m tests.run_all."""

from pathlib import Path
import subprocess
import sys
from time import perf_counter


def main():
    test_directory = Path(__file__).resolve().parent
    root = test_directory.parent
    modules = [f'tests.{path.stem}'
               for path in sorted(test_directory.glob('test_*.py'))]
    if not modules:
        print('No test modules found.', file=sys.stderr)
        return 1

    results = []
    started = perf_counter()
    try:
        for index, module in enumerate(modules, 1):
            print(f'\n[{index}/{len(modules)}] {module}', flush=True)
            module_started = perf_counter()
            try:
                process = subprocess.run(
                    [sys.executable, '-B', '-u', '-m', module], cwd=root,
                )
                code = process.returncode
            except OSError as exc:
                print(f'Could not start {module}: {exc}', file=sys.stderr, flush=True)
                code = 1
            elapsed = perf_counter() - module_started
            results.append((module, code, elapsed))
            status = 'PASS' if code == 0 else f'FAIL (exit {code})'
            print(f'{status}: {module} ({elapsed:.2f}s)', flush=True)
    except KeyboardInterrupt:
        print('\nTest run interrupted.', file=sys.stderr, flush=True)
        return 130

    print('\nTest results', flush=True)
    for module, code, elapsed in results:
        status = 'PASS' if code == 0 else f'FAIL (exit {code})'
        print(f'  {status}: {module} ({elapsed:.2f}s)', flush=True)
    failed = sum(code != 0 for _, code, _ in results)
    print(f'\n{len(results)-failed} passed, {failed} failed; '
          f'total {perf_counter()-started:.2f}s', flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
