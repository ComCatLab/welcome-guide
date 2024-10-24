# Setting Up for Local Development

This page only covers the basic setup to enable local development of the
welcome guide. Building the documentation is covered [elsewhere](./docs.md).

## Installing Prerequisites

1. Install [Python 3.10][python] or greater.

2. (Optional) Install [Hatch][install-hatch]. Handles environment management
   and provides some handy command-line functions to simplify documentation
   building.

3. (Optional) Install [Pandoc][pandoc]. Only required to build the PDF version
   of the guide.

## Getting the Source Files

1. Create your own fork of the [Welcome Guide][welcome-guide] (look for the
   "Fork" button on GitHub).

2. Clone your fork.

    with `ssh` (requires that your `ssh` key is added to your GitHub account):

    ```bash
    git clone git@github.com:yourusername/welcome-guide.git
    ```

    with `https` (requires login to GitHub account):

    ```bash
    git clone https://github.com/yourusername/welcome-guide.git
    ```

3. Initialize the development virtual environment.

    with [Hatch][install-hatch]:

    ```bash
    hatch env create default
    ```

    with `venv` and [pip][pip]:

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    python3 -m pip install '.[test,dev,docs]'
    ```

4. Install the pre-commit hooks.

    ```bash
    pre-commit install-hooks
    ```

[python]: https://www.python.org/downloads/
[install-hatch]: https://hatch.pypa.io/latest/install/
[welcome-guide]: https://github.com/ComCatLab/welcome-guide
[pip]: https://pip.pypa.io/en/stable/
[pandoc]: https://pandoc.org/installing.html
