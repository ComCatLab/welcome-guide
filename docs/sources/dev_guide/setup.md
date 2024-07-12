# Setting Up for Local Development

This page covers the basic setup to be able to locally develop the welcome
guide.

## Installing Prerequisites

1. Install [Python 3.10][python] or greater.

2. (Optional) Install [Hatch][install-hatch]. Optionally handles environment management.

3. (Optional) Install [Pandoc][pandoc]. Only required to build the PDF version of the guide.

## Getting the Source Files

1. Create your own fork of the [Welcome Guide][welcome-guide] (look for the "Fork" button on GitHub).

2. Clone your fork.

    with `ssh` (requires that your `ssh` key is added to your GitHub account):

    ```bash
    git clone git@github.com:yourusername/welcome-guide.git
    ```

    with `https` (requires login to GitHub account):

    ```bash
    git clone https://github.com/yourusername/welcome-guide.git
    ```

3. (Optional) Initialize the development virtual environment.

    with [Hatch][install-hatch]:

    ```bash
    hatch env create default
    ```

    with [pip][pip]:

    ```bash
    python3 -m venv .venv && python3 -m pip install .
    ```

[python]: https://www.python.org/downloads/
[install-hatch]: https://hatch.pypa.io/latest/install/
[welcome-guide]: https://github.com/ComCatLab/welcome-guide
[pip]: https://pip.pypa.io/en/stable/
[pandoc]: https://pandoc.org/installing.html
