# Setting Up the Tests

This page explains the testing framework and explains how to run the tests.

## Testing Framework

This project uses [`pytest`][pytest] to collect and run tests. These tests can
be used to validate that your virtual environment is set up correctly. `pytest`
provides useful features that make it simpler to write tests. Some basic
features are covered below.

### Writing Tests

The following snippet highlights useful principles to follow when writing tests

```python
from pathlib import Path

def test_should_create_existing_temporary_directory(tmp_path: Path) -> None:
    assert tmp_path.exists()
```

Some notes:

- This test check tests that the temporary directory provided exists. If the
  directory does not exist, the assertion will raise an error and `pytest` will
  report that the test fails.
- This test requests the [`tmp_path`][tmp-path] fixture.
- **Each test should be contained in a function**.
- **Tests should be relatively simple** (this is the idea of "unit-testing").
- **Test functions should have descriptive names** of the form
  `test_should_do_something` (where `do_something` is replaced with the
  expected behaviour of the test). It is okay if the name is long or verbose.
  Clarity is the main priority. Someone should be able to get an idea of what
  the function is testing just from reading the name without having to read the
  code. If you're finding, however, that you're not able to simply describe
  what a function should test, that may be a sign that it needs to be broken up
  into smaller tests.
- Strictly speaking, the return annotation for tests should be `None`.
- **The condition for a test should be checked with an `assert` statement.**

### Using `pytest` Fixtures

[`pytest` Fixtures][fixtures] are functions that return a value or perform an action
and can be "requested" by test functions. By "requested", we mean that a test
function can require that the fixture function is executed and, optionally,
that its output be provided to the test function. For example, a fixture may
create an `Atoms` object that is used by a test function to perform a DFT
calculation. Additionally, fixtures may perform tasks like changing the current
working directory, creating files, or requesting other fixtures. `pytest` and
its related plugins provide [built-in fixtures][pytest-fixtures], or you can
[write your own][writing-fixtures]. Fixtures are declared by using the
`pytest.fixture` decorator:

```python
import pytest

from ase.build import molecule

@pytest.fixture(name="co2")
fixture_co2() -> Atoms:
    co2 = molecule("CO2")

    return co2
```

You can then "request" the fixture in a test function like so:

```python
def test_should_create_co2_molecule(co2: Atoms) -> None:
    assert co2.symbols == "CO2"
```

!!! Tip

    If you notice that you are using the same variable in multiple tests, this
    is a great use-case for fixtures! For example:

    ```python
    from ase.build import fcc111

    def test_should_contain_nickel() -> None:
        nickel_unit_cell = fcc111(
            "Ni", size=[1, 2, 1], a=3.55, orthogonal=True
        )
        assert "Ni" in nickel_unit_cell.symbols
    
    def test_should_have_positive_lattice_constant() -> None:
        nickel_unit_cell = fcc111(
            "Ni", size=[1, 2, 1], a=3.55, orthogonal=True
        )
        assert (nickel_unit_cell.cell[0] @ nickel_unit_cell.cell[0])**(1/2) > 0
    ```
    Both test functions request the `nickel_unit_cell` fixture and perform
    separate checks on the returned `Atoms` object. But the same code can be
    simplified as follows:

    ```python
    from ase.build import fcc111

    @pytest.fixture(name="nickel_unit_cell")
    fixture_nickel_unit_cell() -> Atoms:
        nickel_unit_cell = fcc111(
            "Ni", size=[1, 2, 1], a=3.55, orthogonal=True
        )
        return nickel_unit_cell
    
    def test_should_contain_nickel(nickel_unit_cell: Atoms) -> None:
        assert "Ni" in nickel_unit_cell.symbols
    
    def test_should_have_positive_lattice_constant(nickel_unit_cell: Atoms) -> None:
        assert (nickel_unit_cell.cell[0] @ nickel_unit_cell.cell[0])**(1/2) > 0
    ```

    Now, it is more clear what each individual test is testing, and scaling
    the test suite to include more test cases doesn't require repeating the
    same "boilerplate" code.

### Skipping Tests

Tests can be skipped conditionally using [`pytest.mark.skipif`][skip-if].

```python
import sys

import pytest

@pytest.mark.skipif(
    sys.version_info < (3,11),
    reason="requires Python 3.11 or higher"
)
def test_should_be_skipped_if_python_version_less_than_311() -> None:
    pass
```

or at the command line by (de-)selecting tests using [markers][markers].

```shell
pytest -m "not calculator" tests
```

The above snippet deselects all tests marked with the "calculator" marker.

## Running Tests

1. To set up the test framework, ensure that the optional "test" dependencies
   have been installed. Run the following from the root directory of the
   welcome guide to install the required dependencies into the current
   environment.

    ```python
    pip install '.[test]'
    ```

2. Call `pytest`.

    ```shell
    pytest tests
    ```

3. Verify that all tests pass.

    ```shell
    ====================================== test session starts ======================================
    platform darwin -- Python 3.11.9, pytest-8.2.2, pluggy-1.5.0 -- ~/welcome-guide/bin/python
    cachedir: .pytest_cache
    rootdir: ~/welcome-guide
    configfile: pyproject.toml
    plugins: datadir-1.5.0, cov-5.0.0, xdist-3.6.1
    8 workers [5 items]
    scheduling tests via LoadScheduling

    tests/test_gpaw.py::test_gpaw_calculator
    tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer
    tests/test_gpaw.py::test_should_perform_gpaw_calculation
    tests/test_vasp.py::test_vasp_calculator
    [gw3] [ 20%] PASSED tests/test_gpaw.py::test_should_perform_gpaw_calculation
    [gw4] [ 40%] PASSED tests/test_vasp.py::test_vasp_calculator
    [gw1] [ 60%] PASSED tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer
    tests/test_atom_manipulations.py::test_should_create_nickel_slab
    [gw0] [ 80%] PASSED tests/test_atom_manipulations.py::test_should_create_nickel_slab
    [gw2] [100%] PASSED tests/test_gpaw.py::test_gpaw_calculator

    ============================================ PASSES =============================================
    ==================================== short test summary info ====================================
    PASSED tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer
    PASSED tests/test_atom_manipulations.py::test_should_create_nickel_slab
    PASSED tests/test_gpaw.py::test_should_perform_gpaw_calculation
    PASSED tests/test_vasp.py::test_vasp_calculator
    PASSED tests/test_gpaw.py::test_gpaw_calculator
    ====================================== 5 passed in 351.31s ======================================
    ```

    !!! Warning

        Some tests will be skipped if the required dependencies are not installed.
        Others will fail if your environment is not set up for the ASE calculator
        interface.

        ```shell hl_lines="12 14 16 25-27"
        ==================================================== test session starts =====================================================
        platform darwin -- Python 3.11.9, pytest-8.2.2, pluggy-1.5.0 -- ~/welcome-guide/bin/python
        cachedir: .pytest_cache
        rootdir: ~/welcome-guide
        configfile: pyproject.toml
        plugins: datadir-1.5.0, cov-5.0.0, xdist-3.6.1
        8 workers [5 items]     
        scheduling tests via LoadScheduling

        tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer 
        tests/test_vasp.py::test_vasp_calculator 
        [gw4] [ 20%] SKIPPED tests/test_vasp.py::test_vasp_calculator 
        tests/test_gpaw.py::test_gpaw_calculator 
        [gw2] [ 40%] SKIPPED tests/test_gpaw.py::test_gpaw_calculator 
        tests/test_gpaw.py::test_should_perform_gpaw_calculation 
        [gw3] [ 60%] SKIPPED tests/test_gpaw.py::test_should_perform_gpaw_calculation 
        tests/test_atom_manipulations.py::test_should_create_nickel_slab 
        [gw0] [ 80%] PASSED tests/test_atom_manipulations.py::test_should_create_nickel_slab 
        [gw1] [100%] PASSED tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer 

        =========================================================== PASSES ===========================================================
        ================================================== short test summary info ===================================================
        PASSED tests/test_atom_manipulations.py::test_should_create_nickel_slab
        PASSED tests/test_atom_manipulations.py::test_should_create_slab_with_water_layer
        SKIPPED [1] tests/test_vasp.py:14: The environment is not configured for VASP. See the ASE documentation for details.
        SKIPPED [1] tests/test_gpaw.py:13: GPAW is not installed.
        SKIPPED [1] tests/test_gpaw.py:32: GPAW is not installed.
        ================================================ 2 passed, 3 skipped in 2.67s ================================================
        ```

        Follow the instructions in the "short test summary info" to resolve
        the issue.

[tmp-path]: https://docs.pytest.org/en/latest/reference/reference.html#std-fixture-tmp_path
[pytest]: https://docs.pytest.org/en/latest/index.html
[fixtures]: https://docs.pytest.org/en/latest/reference/reference.html#fixtures-api
[pytest-fixtures]: https://docs.pytest.org/en/latest/reference/fixtures.html#reference-fixtures
[writing-fixtures]: https://docs.pytest.org/en/latest/how-to/fixtures.html
[skip-if]: https://docs.pytest.org/en/latest/reference/reference.html#pytest-mark-skipif
[markers]: https://docs.pytest.org/en/latest/example/markers.html#mark-examples
