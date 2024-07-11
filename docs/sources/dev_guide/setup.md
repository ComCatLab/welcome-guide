# Getting Started for Local Development

This page covers the basic setup to be able to locally develop the welcome
guide.

## Installing Prerequisites

1. Install [Python 3.10][python] or greater.

2. Install [Hatch][hatch], the build backend and environment manager for the welcome guide. If
   [`pipx`][pipx] is installed:

   ```bash
   pipx install hatch
   ```

   Otherwise, follow the instructions [here][install-hatch].

## Getting the Source Files

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

## Making Your Changes

1. Change your current working directory to be that of the project root.

    ```bash
    cd welcome-guide
    ```

2. Initialize the development virtual environment.

    with [Hatch][install-hatch]:

    ```bash
    hatch env create default
    ```

    with [pip][pip]:

    ```bash
    python -m venv .venv && python -m pip install .
    ```

3. Create and checkout a branch for your local changes. (Not directly committing changes to the `main` branch
   helps ensure that, if synced with the fork source, everyone's main branch `main` means the same
   thing.)

    ```bash
    git checkout -b name-of-feature-or-change
    ```

4. Activate the development environment.

    If installed via [Hatch][install-hatch]:

    ```bash
    hatch shell
    ```

    If installed via [pip][pip]:

    ```bash
    source .venv/bin/activate
    ```

5. Make your changes (e.g., add/change/remove files).
6. Commit your changes.

    ```bash
    git commit -S -m "I made a change"
    ```

7. Push your changes to **your** remote.

    ```bash
    git push origin name-of-feature-or-change
    ```

8. [Create a pull request][pull-requests] for your change on GitHub.

[python]: https://www.python.org/downloads/
[hatch]: https://hatch.pypa.io/latest/
[install-hatch]: https://hatch.pypa.io/latest/install/
[welcome-guide]: https://github.com/ComCatLab/welcome-guide
[pipx]: https://github.com/pypa/pipx
[pip]: https://pip.pypa.io/en/stable/
[pull-requests]: https://github.com/ComCatLab/welcome-guide/pulls
