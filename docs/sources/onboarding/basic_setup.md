# Basic Setup

## System Requirements

- MacOS 26.0.1 or later

## Prerequisites

These instructions do not assume any familiarity with shell scripting
or Python; however, the following reference pages may be of use:

- [The Z Shell Manual](https://zsh.sourceforge.io/Doc/Release/zsh_toc.html)
- [The Python Tutorial](https://docs.python.org/3/tutorial/index.html)

## Objectives

- to install [miniconda][miniconda]
- to develop a basic understanding of virtual environments
- to create a general purpose virtual environment for use with ComCat Lab
- to understand the `.zshrc` file

## Step-by-Step Instructions

1. **Open the "Terminal" application.**

    (TODO: explain a shell (zsh and bash) and how to navigate in home directory)

2. **Install the XCode Development Tools.**

    ```shell
    xcode-select --install
    ```

    You should see a pop up asking you to download the XCode Development Tools.

    (TODO: add screenshot of pop up)

    Click "Install" and accept the license agreement.

3. **Install [miniconda][miniconda].**

    ```shell
    cd ~/
    mkdir ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    ```

    ??? info "Explanation"

        The first command changes the current working directory to the
        directory given by the argument `~/`. `~/` is a special expression that
        always evaluates to your home directory. (This process is called
        [filename expansion][zsh-filename-expansion].) Grammatically, `cd` is
        commonly used as a verb (e.g., `cd` into a directory). The second
        command, creates a sub-directory of your home directory called
        "miniconda3". That is, thee `mkdir` command "makes directories".
        The third command downloads the miniconda installer (using the `curl`
        command) and saves it to the file `~/miniconda3/miniconda.sh`. The
        fourth command runs the installer, and the final command deletes it.

    (TODO: briefly explain the conda environment/package manager)

4. **Initialize `conda` for shell interaction and activate the base environment.**

    ```shell
    conda init --all
    source ~/miniconda3/bin/activate
    ```

    ??? info "Explanation"

        The first command adds lines to your `.zshrc` file that ensure that your base conda
        environment is activated when you start a new terminal window. The second command
        activates the base environment.

        (TODO: explain VENVs and their purpose)

    The name of the active virtual environment should now be visible in the shell prompt.

    (TODO: add screenshot of terminal prompt with base environment activated)

5. **Create a virtual environment for your work with ComCat Lab.**

    ```shell
    conda create -n comcat python=3.11.13
    conda activate comcat
    ```

    ??? info "Explanation"

        The first command creates a conda environment called "comcat" with version
        3.11.13 of Python. The second command activates the newly created
        environment.

    Please accept the accept the Terms of Service, and type "y" to proceed.

    (TODO: add screenshot of TOS and terminal prompts)

    You should now see the name of your environment in the prefix of your terminal
    prompt.

    (TODO: add screenshot of terminal prompt with `comcat` environment activated)

6. **Install some necessary packages and verify their installation.**

    ```shell
    pip install ase matplotlib ipython ruff pymatgen matplotlib mp-api
    which ase
    ```

    ??? info "Explanation"

        The first command will install several commonly used packages into
        your "comcat" environment. The second command should result in the
        path to the ASE executable being printed to your terminal, if successful.
        (TODO: explain pip vs. conda)

    You should also confirm that the other packages are "importable" by
    launching a Python interpreter:

    ```shell
    python
    ```

    (TODO: add screenshot of Python interpreter)

    ```python
    import ase, matplotlib, ruff, pymatgen, matplotlib, mp-api
    ```

    If this command fails, then something went wrong with the `pip` command
    above. Exit the Python interpreter by typing:

    ```python
    exit()
    ```

    !!! Tip

        Check out the [software pages](../software_pages.md) for information about
        each package.

7. **Define an alias for launching the ASE GUI.**

    Using a text editor of your choice (e.g., `vim`, `nano`, or `emacs`), open your
    `.zshrc` file for editing. Your `.zshrc` file is run every time a new interactive
    non-login shell is launched. (Generally, every time you open a terminal
    window, you will be launching an interactive non-login shell.)

    ```shell
    vim ~/.zshrc
    ```

    !!! Tip

        Type `vimtutor` into your terminal to launch the `vim` tutorial!

    or equivalently,

    ```shell
    nano ~/.zshrc
    ```

    Add the following line to your `.zshrc` file:

    ```shell
    alias ag='ase gui'
    ```

    ??? info "Explanation"

        Adding this line defines an *alias* so that you can launch the ASE GUI
        with the command `ag` instead of the full command `ase gui` from the
        command line. (`alias` is a [shell builtin command][alias].)

    Close the file (by typing `:wq` in `vim` or typing `ctrl+x` in `nano`), and
    then run the following from the command-line:

    ```shell
    source .zshrc
    ```

    ??? info "Explanation"

        In order for any changes in your `.zshrc` to take effect for the current
        terminal session, you **must** "source" the file. ("Sourcing" a file
        executes the commands in the file. Alternatively, opening
        a new terminal session will automatically source the file.)

[miniconda]: https://www.anaconda.com/docs/getting-started/miniconda/main
[alias]: https://www.gnu.org/software/bash/manual/bash.html#Aliases
[zsh-filename-expansion]: https://zsh.sourceforge.io/Doc/Release/Expansion.html#Filename-Expansion
