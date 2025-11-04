# Local Setup

This page contains walkthroughs for setting up your local machine.

## Basic Setup

### System Requirements

- MacOS 26.0.1 or later

### Objectives

- to install [miniconda][miniconda]
- to develop a basic understanding of virtual environments
- to create a general purpose virtual environment for use with ComCat Lab
- to understand the `.zshrc` file

### Step-by-Step Instructions

These instructions do not assume any familiarity with shell scripting
or Python; however, the following reference pages may be of use:

- [The Z Shell Manual](https://zsh.sourceforge.io/Doc/Release/zsh_toc.html)
- [The Python Tutorial](https://docs.python.org/3/tutorial/index.html)

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

## Software Development Setup

### System Requirements

- MacOS 26.0.1 or later
- A [valid CCDB account](../tutorials/ccdb.md)
- A valid [GitHub][github] account
- [VSCode 1.105.1 or later](https://code.visualstudio.com)

### Objectives

- to configure your Git profile
- to gain familiarity with SSH
- to create an SSH key
- to create your first repository
- to gain familiarity with Git, GitHub, and version control
- to learn what an IDE is

### Step-by-Step Instructions

If you are not familiar with [Git][git] and version control, it is highly recommended
that you read the following short pages:

- [About Version Control](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control)
- [What is Git?](https://git-scm.com/book/en/v2/Getting-Started-What-is-Git%3F)

1. **Open the "Terminal" application.**

2. **Configure your Git profile.**

    ```shell
    git config --global user.name <name>
    git config --global user.email <email>
    git config --global core.editor <editor>
    git config --global init.defaultBranch main
    ```

    ??? info "Explanation"

        These four commands set your Git username, email, preferred text
        editor, and the name of the default branch (created when you initialize
        a repository). The `--global` option is used to ensure that the
        configuration holds for all repositories on your current machine, unless
        otherwise stated. The `--local` option can be used to make repository-
        specific configurations.

    Your name and email will be used to attribute to you changes that you make to
    git repositories. The editor you specify will be launched by git for editing
    commits. For more tips on how to configure git, check out
    [this page][git-config].

    <a name="create-ssh-key"></a>

3. **Create an SSH key.**

    SSH is a network protocol that allows for secure connections between computers.
    It is the primary way that you will connect to Digital Research Alliance of
    Canada clusters and is also a preferred method of interacting with remote git
    repositories. To use SSH, you must create an SSH key, which can be done using
    the following command:

    ```shell
    ssh-keygen -t ed25519 -C <your_email@here.com> -f ~/.ssh/id_ed25519
    ```

    ??? info "Explanation"

        (TODO: explain command and options)

    Enter a passphrase for your SSH key when prompted. You will later need the
    contents of your public key file. To display the contents of the file, type
    the following command into your terminal and copy the output:

    ```shell
    cat ~/.ssh/id_ed25519
    ```

4. **Add your SSH key to your GitHub profile.**

    [Follow the instructions from GitHub][github-ssh].

    ??? info "Explanation"

        This will enable you to securely connect to the GitHub servers and
        download repositories to which you have access.

5. **Add your SSH key to your CCDB account.**

    Log in via [https://ccdb.alliancecan.ca](https://ccdb.alliancecan.ca).

    Then, under the "My Account" tab, navigate to the "SSH Keys" link.

    Copy the contents of your public SSH key file.

    Give your key a description (e.g., ComCatLab iMac)

    Press "Add Key".

6. **Create your first Git repository.**

    Open Terminal and navigate to the directory in which you would
    like to create your repository.

    ```shell
    mkdir -p ~/Repositories/comcat-tutorial
    cd ~/Repositories/comcat-tutorial
    git init
    ```

    ??? info "Explanation"

        First, we create the directory that will house our repository and
        `cd` into it. Then, by calling `git init`, we initialize the repository
        with `git`.

7. **Create a `.gitignore` file.**

    ```shell
    vim .gitignore
    ```

    ??? info "Explanation"

        This file defines filenames and patterns that Git should ignore. This is
        especially helpful if there are very large files that you don't want to
        track with Git or if there are files with proprietary/sensitive
        information that you don't want to accidentally check into version control.

    Add the following lines to the `.gitignore` file and close the file:

    ```shell
    .vscode
    __pycache__
    ```

    ??? info "Explanation"

        We ignore files named `.vscode` because this is the name of a
        configuration directory for a software called VSCode. Settings
        defined in files within this directory are user-specific and so
        should not be tracked by Git. `__pycache__` is a generic folder
        name used by Python to store byte-code generated from running
        scripts in order to speed up successive executions of the same
        code. These folders are generated every time you run Python code
        based on the actual Python code in files. They are system-specific
        and thus should not be tracked with Git.

8. **Create your first commit**.

    To record your changes with Git, you must first add the edited files
    to the index using `git add`.

    ```shell
    git add .gitignore
    ```

    !!! Tip

        You can add all files in the current directory to the index using
        `git add .`

    Saving (or committing) your changes with Git is accomplished by creating a
    "commit". This creates a snapshot that can be used to rollback changes
    or apply changes to another branch.

    ```shell
    git commit
    ```

    This will launch the editor that you specified in Step 2 for you to write
    a commit message. **A good commit message features a header line
    that summarizes the changes in 80 characters or less.** A more comprehensive
    description/explanation can be included below (separated from the header
    by a blank line). For example, a reasonable message for the current commit
    could be:

    ```shell
    Created .gitignore file

    The .vscode configuration directory and Python byte-code files are now
    ignored.
    ```

     Once you have written a satisfactory message,
    saving and closing the file will automatically complete the commit.

    !!! Tip

        - The `-m` option can be used to specify the commit message from the
          command-line.
        - The `--amend` option can be used to correct the commit message of the
          previous commit.
        - Saving and closing an empty file will abort the commit.

        Check out [this article][git-best-practices] for guidance on Git best
        practices.

9. **Create a GitHub repository.**

    Before you can push your changes to GitHub, you will first need to create a
    repository on GitHub. First, log in to your GitHub account.

    Click on your profile icon and navigate to the "Repositories" tab.

    Click on "New" to create a new repository.

    Name and describe your repository.

    Configure your repository. For now, set the visibility to "Private",
    do not add a .gitignore or license (we will do this ourselves), but
    check "Yes" for "Add a README".

    Click "Create Repository".

10. **Make changes to your repository on GitHub.**

    (TODO: Write README description)

    ??? info "Explanation"

        (TODO: allude to other workflows for making changes (issues, PRs))

11. **Sync changes from a remote repository.**

    Note that since we have files in the GitHub repository that are not
    in your local version of the repository (and vice-versa), the two copies
    of the repository are out of sync. To sync them, we must first add the
    GitHub repository as a "remote". From a directory in your repository, run:

    ```shell
    git remote add origin git@gitlab.com:<your-github-username>/<gh-repo>
    ```

    ??? info "Explanation"

        (TODO: explain command and link relevant git book section for
        distributed Git on remotes)

    We can then "fetch" the changes from the remote using the following:

    ```shell
    git fetch
    ```

    ??? info "Explanation"

        (TODO: explain the concept of fetching)

    We can merge the changes from the remote into our local branch with
    the command

    ```shell
    git merge origin/main
    ```

    ??? info "Explanation"

        (TODO: explain the concept of fetching)

    !!! Tip

        The process of fetching and merging can be automated with the
        command `git pull`. This command will automatically fetch the
        changes from the given remote and merge them into the current
        branch. Note that an upstream branch must be set prior using
        the `--set-upstream` (or `-u`) option.

12. **Push your changes to GitHub.**

    Then, we can push our local changes to GitHub.

    ```shell
    git push --set-upstream origin git@gitlab.com:<your-github-username>/<gh-repo>
    ```

    ??? info "Explanation"

        (TODO: explain command and how both changes are visible on GH)

    !!! Tip

        For a more in-depth explanation of the concepts of committing,
        branching, and merging checkout the [Git Basics][git-basics] and
        [Git Branching][git-branching] sections in the Git book.

13. **Configure VSCode.**

    (TODO: install extensions, create venv, requirements.txt, pre-commit.yaml,
    scripts, .gitignore, connect GH account)

[miniconda]: https://www.anaconda.com/docs/getting-started/miniconda/main
[alias]: https://www.gnu.org/software/bash/manual/bash.html#Aliases
[github]: http://github.com
[git]: https://git-scm.com/book/en/v2
[git-vcs]: https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control
[git-config]: https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration
[github-ssh]: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
[zsh-filename-expansion]: https://zsh.sourceforge.io/Doc/Release/Expansion.html#Filename-Expansion
[git-best-practices]: https://about.gitlab.com/topics/version-control/version-control-best-practices/
[git-basics]: https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository
[git-branching]: https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell
