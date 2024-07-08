# Getting Started for Local Development

This page covers the basic setup to be able to locally develop the welcome
guide.

## Installing Pre-requisites

1. Install [Python 3.10][python] or greater.

2. Install [Hatch][install-hatch] the build backend and environment manager for the welcome guide. If
   [`pipx`][pipx] is installed:

   ```bash
   pipx install hatch
   ```

## Setup

1. Create your own fork of the [Welcome Guide][welcome-guide] (look for the "Fork" button on GitHub).

2. Clone your fork [Welcome Guide repository][welcome-guide].

   with `ssh` (requires that your `ssh` key is added to your GitHub account):

    ```bash
    git clone git@github.com:yourusername/welcome-guide.git
    ```

   with `https` (requires login to GitHub account):

    ```bash
    git clone https://github.com/yourusername/welcome-guide.git
    ```

3. Create and checkout a branch for your local changes. (Not directly committing changes to the `main` branch
   helps ensure that, if synced with the fork source, everyone's main branch `main` means the same
   thing.)

   ```bash
   git checkout -b name-of-feature-or-change
   ```

4. Install the development virtual environment.

   with [Hatch][install-hatch]:

   ```bash
   cd welcome-guide && hatch env create default
   ```

   with [pip][pip]:

   ```bash
   cd welcome-guide && python -m venv .venv && python -m pip install .
   ```

5. Activate the development environment.

   If installed via [Hatch][install-hatch]:

   ```bash
   hatch shell
   ```

   If installed via [pip][pip]:

   ```bash
   source .venv/bin/activate
   ```

[python]: https://www.python.org/downloads/
[install-hatch]: https://hatch.pypa.io/latest/install/
[welcome-guide]: https://github.com/ComCatLab/welcome-guide
[pipx]: https://github.com/pypa/pipx
[pip]: https://pip.pypa.io/en/stable/
