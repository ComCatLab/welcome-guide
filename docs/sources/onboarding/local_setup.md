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

2. Install the XCode Development Tools

    ```shell
    xcode-select --install
    ```

3. Install [miniconda][miniconda].

    ```shell
    cd ~/
    mkdir -p ~/miniconda3
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm ~/miniconda3/miniconda.sh
    ```

4. Initialize conda for shell interaction and activate the miniconda base environment.

    ```shell
    conda init –all
    source ~/miniconda3/bin/activate
    ```

    (TODO: briefly explain virtual environments and the conda environment/package manager)
    (TODO: add screenshot of terminal prompt with base environment activated)

5. Create an environment for your work with ComCat Lab.

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

6. Install some necessary packages.

    ```shell
    conda install ase matplotlib ipython ruff pymatgen matplotlib mp-api
    conda activate comcat
    ```

    !!! Tip

        Check out the software page for information about each package.

7. Edit your `.zshrc` file to configure your shell environment using `vim`.

    ```shell
    cd ~/
    vim .zshrc
    ```

    Your `.zshrc` file is a start-up file-it is run every time a new interactive
    non-login shell is launched. Generally, every time you open a terminal
    window, you will be launching an interactive non-login shell.

    !!! Tip

        Type `vimtutor` into your terminal to launch the vim tutorial!

8. Add the following lines to your `.zshrc` file:

    ```shell
    alias ag='ase gui'
    ```

    Adding this line defines an *alias* so that you can launch the ASE GUI
    with the command `ag` instead of the full command `ase gui` from the
    command line. (`alias` is a [shell builtin command][alias].)

9. Close the file by typing `:wq`, and then run the following from the command-
   line:

    ```shell
    source .zshrc
    ```

    !!! Tip

        In order for any changes in your `.zshrc` to take effect for the current
        terminal session, you **must** source the file. Alternatively, opening
        a new terminal session will automatically source the file.

### Remote Development Setup

If you are not familiar with `git` and version control, it is highly recommended
that you read the following short pages:

- [About Version Control](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control)
- [What is Git?](https://git-scm.com/book/en/v2/Getting-Started-What-is-Git%3F)

1. Open the "Terminal" application.

2. Configure your [git][git] profile with your name, your SFU email, and your
   preferred text editor.

    ```shell
    git config --global user.name <name>
    git config --global user.email <email>
    git config --global core.editor <editor>
    ```

    Your name and email will be used to attribute to you changes that you make to
    git repositories. The editor you specify will be launched by git for editing
    commits. For more tips on how to configure git, check out
    [this page][git-config].

    <a name="create-ssh-key"></a>

3. Next, you will create an SSH key.SSH is a network protocol that allows for secure
   connections between computers. It is the primary way that you will connect to
   Digital Research Alliance of Canada clusters and is also a preferred method
   of interacting with remote git repositories.

    ```shell
    ssh-keygen -t ed25519 -C <your_email@here.com>
    ```

    You will be prompted to enter the file in which to save the key. Press
    return to accept the default location. Next, copy the contents of your
    public key file. To display the contents of the file, type:

    ```shell
    cat /your/ssh_key/location
    ```

4. Add your SSH key to your GitHub profile by following the instructions
   [from GitHub][github-ssh].

5. Add your SSH key to your CCDB account by first logging in via
   [https://ccdb.alliancecan.ca](https://ccdb.alliancecan.ca).

   Then, under the "My Account" tab, navigate to the "SSH Keys" link.

   Copy the contents of your public SSH key file.

   Give your key a description (e.g., ComCatLab iMac)

   Press "Add Key".

6. Create your first git repository by opening Terminal and navigating to
   the directory in which you would like to create your repository.

   ```shell
   mkdir -p ~/Repositories/comcat-tutorial
   cd ~/Repositories/comcat-tutorial
   git init
   ```

   First, we create the directory that will house our repository, and then,
   by calling `git init`, we initialize the repository with `git`.

[miniconda]: https://www.anaconda.com/docs/getting-started/miniconda/main
[alias]: https://www.gnu.org/software/bash/manual/bash.html#Aliases
[git]: https://git-scm.com/book/en/v2
[git-vcs]: https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control
[git-config]: https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration
[github-ssh]: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
