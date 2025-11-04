# Local-Setup

This page will walk you through setting up your local machine.

## Pre-Requisites

- MacOS 26.0.1 or later
- A [valid CCDB account](/docs/sources/tutorials/ccdb.md)

## Objectives

- to set up your workstation with basic research software
- to develop a basic understanding of virtual environments
- to understand the `.zshrc` file

## Instructions

### Basic Setup

1. Open the "Terminal" application.

2. Install [miniconda][miniconda].

    ```shell
    cd ~/
    mkdir -p ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    ```

3. Initialize conda for shell interaction and activate the miniconda base environment.

    ```shell
    conda init –all
    source ~/miniconda3/bin/activate
    ```

    (TODO: briefly explain virtual environments and the conda environment/package manager)
    (TODO: add screenshot of terminal prompt with base environment activated)

4. Create an environment for your work with ComCat Lab.

    ```shell
    conda create –n comcat python=3.11.13
    conda activate comcat
    ```

    You should now see the name of your environment in the prefix of your terminal
    prompt.

    (TODO: add screenshot of terminal prompt with `comcat` environment activated)

    !!! Tip

        `3.11.9` specifies the version of Python that you would like to use to create
        the virtual environment, and `comcat` can be replaced to rename your environment.

5. Install some necessary packages.

    ```shell
    conda install ase matplotlib ipython ruff pymatgen matplotlib mp-api
    conda activate comcat
    ```

    !!! Tip

        Check out the software page for information about each package.

6. Edit your `.zshrc` file to configure your shell environment using `vim`.

    ```shell
    cd ~/
    vim .zshrc
    ```

    Your `.zshrc` file is a start-up file-it is run every time a new interactive
    non-login shell is launched. Generally, every time you open a terminal
    window, you will be launching an interactive non-login shell.

    !!! Tip

        Type `vimtutor` into your terminal to launch the vim tutorial!

7. Add the following lines to your `.zshrc` file:

    ```shell
    alias ag='ase gui'
    ```

    Adding this line defines an *alias* so that you can launch the ASE GUI
    with the command `ag` instead of the full command `ase gui` from the
    command line. (`alias` is a [shell builtin command][alias].)

8. Close the file by typing `:wq`, and then run the following from the command-
   line:

    ```shell
    source .zshrc
    ```

    !!! Tip

        In order for any changes in your `.zshrc` to take effect for the current
        terminal session, you **must** source the file. Alternatively, opening
        a new terminal session will automatically source the file.

### Remote Development Setup

[miniconda]: https://www.anaconda.com/docs/getting-started/miniconda/main
[alias]: https://www.gnu.org/software/bash/manual/bash.html#Aliases
