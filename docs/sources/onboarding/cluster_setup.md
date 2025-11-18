# Cluster Setup

## System Requirements

- MacOS 26.0.1 or later
- A [valid CCDB account](../tutorials/ccdb.md)
- [VSCode 1.105.1 or later](https://code.visualstudio.com)

## Prerequisites

Although this tutorial does not require a deep understanding of the following
concepts, the following documentation pages may prove useful:

- [`scp`](https://docs.alliancecan.ca/wiki/Transferring_data#SCP)
- [`ssh`](https://docs.alliancecan.ca/wiki/SSH)
- [Linux introduction](https://docs.alliancecan.ca/wiki/Linux_introduction)

## Objectives

- to copy files with `scp`
- to login into the login node of Digital Research Alliance (DRA) of Canada
  clusters
- to configure your DRA cluster accounts for running calculations

## Step-by-Step Instructions

This tutorial uses the Python package, `cluster-setup` to template
and configure your DRA cluster directories. By defining a static configuration
file, you can automate the setup of your directories and software creation
across multiple DRA clusters. This following steps outline how to perform
this setup on Fir, but setup on any one of the other clusters is analogous.

1. **Download the following files for setting up your cluster account.**

    - Software scripts: With regards to `cluster-setup`, "software scripts" are
      executable files that install and configure software. They are specified
      with the `--software-script` CLI option or `software_scripts` key in the
      `cluster-setup` configuration file.
        - `install_aliases.bash`: This script will install a file that defines
          aliases that can be encapsulated in a module.
        - `install_shell_customization.bash`: This script
        - `configure_software.bash`: This script configures [ASE][ase],
          [autojob][autojob], and [ccu][ccu] for use by creating suitable
          configuration files
        - `configure_vasp.bash`: This script configures VASP to be used by ASE
          by copying various support files (pseudopotentials, vDW-DF kernel, and a
          Python script used to call VASP) to a directory.
    - Modulefile templates: Jinja2 Lua modulefile templates to be used to define
      custom modules on DRA clusters.
        - `aliases/0.0.1.lua.j2`: A modulefile template corresponding to
          `install_aliases.bash`
        - `shell-customization/0.0.1.lua.j2`: A modulefile template corresponding to
          `install_shell_customization.bash`

    !!! note

        Reach out to a group member for instructions of how to get set up with
        the VASP files.

2. **Copy the files to a DRA cluster.**

    First, collect the files into a directory.

    ```shell
    mkdir cluster-setup-files
    cp FILE_1 FILE_2 FILE_3 ... cluster-setup-files
    ```

    Then, copy the files to a DRA cluster using `scp`.

    ```shell
    scp -r cluster-setup-files <dra_username>@fir.alliancecan.ca
    ```

3. **`ssh` into the login node of the cluster.**

    To `ssh` into Fir, type:

      ```shell
      ssh -Y <dra-username>@fir.alliancecan.ca
      ```

    with `<dra-username>` replaced by your Digital Research Alliance username,
    and authenticate with Two-Factor authentication.

    ??? tip

        It is conveninent to define aliases for `ssh`ing into each of the
        DRA clusters in your `~/.zshrc` file. To do so, first define a variable
        `DRA_USER` like so:

        ```shell
        DRA_USER=<dra-username>
        ```

        Then copy and execute the following command into your terminal:

        ```shell
        for host in fir killarney tamia vulcan nibi narval rorqual trillium; do
        echo "alias $host=ssh -Y $DRA_USER$@$host.alliancecan.ca" >> ~/.zshrc
        done
        source ~/.zshrc
        ```

        `ssh`ing into any cluster can now be accomplished by simply typing its
        name into your terminal:

        ```shell
        fir
        ```

    You should see your newly copied folder in your home directory:

    ```shell
    ls cluster-setup-files
    ```

4. **Create a virtual environment in which to install the `cluster-setup` package.**

    ```shell
    python -m venv cluster-setup-files/.venv && source cluster-setup-files/.venv/bin/activate
    pip install cluster-setup[test]
    ```

    ??? info "Explanation"

        The first command creates a Python virtual environment at `.venv` in the
        current directory and activates it. The final command installs the
        `cluster-setup` package and its `test` extra. Installation
        of the `test` extra enables one to run the test suite prior in
        executing `cluster-setup`.

5. **Run the `cluster-setup` test suite.**

    ```shell
    cluster-setup --test
    ```

    !!! tip

        On Rorqual, it is a known bug that running the tests with the above
        command results in an error related to thread creation. To remedy
        this error, run the test suite with the additional argument to the
        `--test` option.

        ```shell
        cluster-setup --test '-n 64'
        ```

        This limits the number of cores used to 64.

    A test report will be written to the current directory. If the test
    suite runs successfully, you should receive a notification that all tests
    have passed. If not, please open a GitHub issue and attach the test
    report.

6. **Create a configuration file.**

    `cluster-setup` provides a utility for generating a stub configuration file.
    To generate a minimal configuration file that can be filled in, run:

    ```shell
    cluster-setup --config-gen
    ```

    At the very least, edit the following values:

    - `python_packages`: A list of strings indicating Python packages to install
      in your home environment. The same syntax used by `pip install` is
      supported here. For example:

      ```toml
      ...
      python_packages = [
          "ase>=3.25.0",
          "pymatgen",
          "FireWorks",
          "maggma",
      ]
      ...
      ```

    - `git_user_name`: The name with which to sign-off on Git commits
    - `git_email`: The email to associate with your Git commits
    - `git_editor`: The editor to launch when writing commit messages
    - `software_scripts`: A list of software script specifications that is used
      to install software. A software script spec is defined by
      `SCRIPT[:[TEMPLATE]:[MODULE]:[VERSION]:[ARGS]]`. For a detailed
      description, read the `cluster-setup` help text (type `cluster-setup -h`).
      For example, to use the `configure_software.bash` and `configure_vasp.bash`
      software scripts, add something like the following to your configuration
      file:

      ```toml
      ...
      software_scripts = [
          "path/to/configure_software.bash::::{support_file_home} ~/.config",
          "path/to/configure_vasp.bash::::{support_file_home}",
      ]
      ...
      ```

7. **Execute `cluster-setup`.**

    ```shell
    cluster-setup --config-file=config.toml
    ```

    !!! note

        The above command assumes that you have saved your configuration file
        in the current working directory. `config.toml` should be replaced
        by the path to the configuration file.

    This command can take up to five minutes to execute on some clusters. Once
    complete, you should receive a notification. Upon success, deactivate the
    virtual environment

    ```shell
    deactivate
    ```

8. **Verify the setup.**

    Source your login file

    ```shell
    source ~/.bashrc
    ```

    Try to activate your Python environment:

    ```shell
    activate_env
    which python
    ```

    The path that is printed to your terminal should be in a subdirectory
    of your `~/software` directory (or whatever value you entered for the
    `--software-home` option or `software_home` configuration value).

9. **Clean up.**

    The `cluster-setup-files/` directory can now be safely deleted.

    ```shell
    rm -rf cluster-setup-files/
    ```

10. **Rinse and repeat.**

    Log out of the cluster

    ```shell
    exit
    ```

    and repeat steps 2-9 on all clusters on which you would like to run
    calculations.

[ase]: http://ase-lib.org
[autojob]: https://python-autojob.readthedocs.io/en/development/
[ccu]: https://python-comp-chem-utils.readthedocs.io/en/development/
