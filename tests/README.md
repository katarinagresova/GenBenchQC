# Tests

For testing, you need to install `pytest` package.

To run a specific test

```bash
    pytest -v ./tests/test_specific_file.py
```

To run all tests

```bash
    pytest -v tests/
```

## How to write tests

The structure of the test folder mirrors the package structure (note the prefix `test_`). Try to cover both the typical and the corner cases.