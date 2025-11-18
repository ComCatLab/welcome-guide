# Software Development Setup

## System Requirements

- MacOS 26.0.1 or later
- A [valid CCDB account](../tutorials/ccdb.md)
- A valid [GitHub][github] account
- [VSCode 1.105.1 or later](https://code.visualstudio.com)

## Prerequisites

If you are not familiar with [Git][git] and version control, it is highly recommended
that you read the following short pages:

- [About Version Control](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control)
- [What is Git?](https://git-scm.com/book/en/v2/Getting-Started-What-is-Git%3F)

## Objectives

- to configure your Git profile
- to gain familiarity with SSH
- to create an SSH key
- to create your first repository
- to gain familiarity with Git, GitHub, and version control
- to learn what an IDE is

## Step-by-Step Instructions

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
    ssh-keygen -t ed25519 -C <your_email@here.com>
    ```

    Accept the default file location for the private key (usually,
    `~/.ssh/id_25519`) and enter a passphrase for your SSH key when prompted.
    You will later need the contents of your **public** key file. To display
    the contents of the file, type the following command into your terminal
    and copy the output:

    ```shell
    cat ~/.ssh/id_ed25519.pub
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
    mkdir -p ~/Repositories/comcat-project
    cd ~/Repositories/comcat-project
    git init
    ```

    ??? info "Explanation"

        First, we create the directory that will house our repository and
        `cd` into it. Then, by calling `git init`, we initialize the repository
        with `git`. `comcat-project` will be the name of the repository; you
        can replace this name as you wish.

7. **Create a `.gitignore` file.**

    ```shell
    nano .gitignore
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

8. **Create a README file.**

    ```shell
    nano README.md
    ```

    ??? info "Explanation"

        A README file explains the purpose of a directory and the files
        contained therein. It is best practice to create a new README file
        for every project. Doing so ensures that the person returning to your
        project folder (yourself included!) can understand the purpose of the
        calculations therein. Update this file as your project progresses.

    In the file, write a description for your project:

    ```markdown
    # TITLE GOES HERE

    Description goes here...
    ```

    ??? info "Note"

        The file extension `.md` indicates that this file should be interpreted
        as Markdown. Markdown is a markup language that allows you
        add special formatting to otherwise plaintext text documents. Note that
        we have preceded our title with `#`; this syntax is used to denote
        headings. Subsequent subheadings can be created using increasing
        numbers of `#` signs. [Check out this webpage][markdown] for a quick
        primer on Markdown syntax.

9.  **Create your first commit**.

    To record your changes with Git, you must first [start tracking][tracking-files]
    the new files using `git add`.

    ```shell
    git add .gitignore README.md
    ```

    !!! Tip

        You can start tracking all files in the current directory using `git add .`

    Saving (or committing) your changes with Git is accomplished by creating a
    "commit". This creates a snapshot that can be used to rollback changes
    or apply changes to another branch.

    ```shell
    git commit
    ```

    This will launch the editor that you specified in Step 2 for you to write
    a commit message. **A good commit message features a header line
    that summarizes the changes in 50 characters or less.** A more comprehensive
    description/explanation can be included below (separated from the header
    by a blank line). For example, a reasonable message for the current commit
    could be:

    ```shell
    Created .gitignore and README

    The .vscode configuration directory and Python byte-code files are now
    ignored.
    ```

    Once you have written a satisfactory message, saving and closing the file will
    automatically complete the commit.

    !!! Tip

        - The `-m` option can be used to specify the commit message from the
          command-line.
        - The `--amend` option can be used to correct the commit message of the
          previous commit.
        - Saving and closing an empty file will abort the commit.

        Check out [this article][git-best-practices] for guidance on Git best
        practices.

10. **Create a GitHub repository.**

    In order to be able to access your repository from other computers, you
    will need a remote version control system-a central, remote server to store
    your project files and history. This role will be fulfilled by GitHub.
    However, before you can push your changes to GitHub, you will first need to
    create a repository on GitHub. First, log in to your GitHub account.

    Click on your profile icon and navigate to the "Repositories" tab.

    Click on "New" to create a new repository.

    Name and describe your repository. The name should match with the name
    you chose when you first created your repository (e.g., `comcat-project`).

    Configure your repository. For now, set the visibility to "Private",
    do not add a .gitignore, README, or license (we will do this ourselves).

    Click "Create Repository".

    ??? Tip

        You can also edit files within you repository from the GitHub
        interface. For example, try editing your README from your repository
        home.

11. **Sync changes from a remote repository.**

    In general, whenever you restart work on a local version of your repository,
    it is a good idea to ensure that your version of the files in your
    repository (e.g., on your computer) match those in the remote repository
    (on GitHub). This is accomplished by *syncing with the remote*.

    ??? Note

        For example, if you edited your README on GitHub, then the versions
        on GitHub and your computer will not match.

    To sync them, we must first add the GitHub repository as a [remote][remote-repo].
    From a directory in your repository, run:

    ```shell
    git remote add origin git@github.com:<your-github-username>/<gh-repo>.git
    ```

    ??? info "Explanation"

        This command adds your GitHub repository as a remote repository with
        the name `origin`. (This is a common naming scheme used to refer to
        the special remote repository corresponding to one's own version of a
        repository.) `upstream` is another common special name; it is used to
        refer to the original version of a [forked][github-forked] repository.

    !!! note "Important"

        Don't forget to replace `<your-github-username>` and `<gh-repo>` with your
        GitHub username and the repository name, respectively.

    We can then "fetch" the changes from the remote using the following:

    ```shell
    git fetch
    ```

    ??? info "Explanation"

        "Fetching" refers to the act of downloading changes that have been
        made to remote branches since the last "fetch". Note that this does not
        change your version of the branch-it only makes the changes available
        in your local copy of the remote branch. For example, if there have been
        changes to the `main` branch of your `origin` remote, then the changes
        will be available in the `origin/main` branch. But your `main` branch
        will remain untouched until your merge these changes.

    We can merge the changes from the remote into our local branch with
    the command

    ```shell
    git merge origin/main
    ```

    ??? info "Explanation"

        "Merging" changes refers to the act of incorporating the changes from
        one branch into another branch. [Check out this short article][git-merge]
        for a more detailed explanation of merging.

    !!! Tip

        The process of fetching and merging can be automated with the
        command `git pull`. This command will automatically fetch the
        changes from the given remote and merge them into the current
        branch. Note that an upstream branch must be set prior using
        the `--set-upstream` (or `-u`) option of `git push` (see below).

12. **Push your changes to GitHub.**

    Now that we have synced our changes, we can be sure that our version of the
    code only features updates to the repository, and we can push our local
    changes to GitHub.

    ```shell
    git push --set-upstream origin git@github.com:<your-github-username>/<gh-repo>.git
    ```

    ??? info "Explanation"

        The `git push` command is used to send your changes to a remote. Here,
        the `--set-upstream` option is used to specify to which remote the
        changes will be sent. ("upstream" here is not meant to be confused with
        the previous naming convention for remotes.) This option also sets the
        default remote to be used for `git push` commands and can thus be
        omitted for subsequent `git push`es to this remote. For example, pushing
        changes to your `origin` remote can be accomplished by simply typing:

            ```shell
            git push
            ```

    !!! Tip

        For a more in-depth explanation of the concepts of committing,
        branching, and merging checkout the [Git Basics][git-basics] and
        [Git Branching][git-branching] sections in the Git book.

13. **Configure VSCode.**

    VSCode is an integrated development environment (IDE), a software that
    provides features to support software development. Many of the tasks
    completed herein (e.g., repository creation, file editing, syncing with
    remotes) can be accomplished from the VSCode application. VSCode also
    offers several very useful features such as automatic code refactoring,
    advanced search and replace across files, automatic code formatting,
    integration with code linters/formatters, automatic detection of
    virtual environments, and previews of markdown documents.

    To get started, [download and install VSCode][install-vscode],
    [connect your GitHub account][connect-github-vscode], open your newly
    created repository, and install some extensions
    to help the development process. Recommended extensions are listed
    in the [ComCat Lab VSCode Setup Extension Pack][comcat-vscode-extensions].
    At the very least, it is highly recommended that you install the following:

    - [autoDocstring][autodocstring]: Automatically generate Python docstrings
      after typing `"""` and then pressing `<tab>`.
    - [Ruff][ruff-extension]: format and check your Python code for errors
    - [Python][python-extension]: Python language support and debugger
      - This extension enables VSCode to detect syntax errors in your code
        without having to explicitly run it and also comes with a debugger
        that allows you to execute your code step-by-step with access to the
        runtime values of variables in order to identify bugs

    (TODO: install extensions, create venv, requirements.txt, pre-commit.yaml,
    scripts, .gitignore, connect GH account)

[github]: http://github.com
[git]: https://git-scm.com/book/en/v2
[git-vcs]: https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control
[git-config]: https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration
[github-ssh]: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
[git-best-practices]: https://about.gitlab.com/topics/version-control/version-control-best-practices/
[git-basics]: https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository
[git-branching]: https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell
[tracking-files]: https://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository#:~:text=tracking%20the%20file.-,Tracking%20New%20Files,-In%20order%20to
[markdown]: https://www.markdownguide.org/basic-syntax/
[remote-repo]: https://git-scm.com/book/en/v2/Git-Branching-Remote-Branches
[git-merge]: https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging
[github-forked]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo
[install-vscode]: https://code.visualstudio.com/download
[connect-github-vscode]: https://code.visualstudio.com/docs/sourcecontrol/github
[autodocstring]: https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring
[ruff-extension]: https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff
[Python-extension]: https://marketplace.visualstudio.com/items?itemName=ms-python.python
