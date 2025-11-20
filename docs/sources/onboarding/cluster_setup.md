# Cluster Setup

## System Requirements

- MacOS 26.0.1 or later
- A [valid CCDB account](../tutorials/ccdb.md)

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

This tutorial uses the Python package, [`cluster-setup`][cluster-setup] to
template and configure your DRA cluster directories. By defining a static
configuration file, you can automate the setup of your directories and software
creation across multiple DRA clusters. The following steps outline how to
perform this setup on Fir, but setup on any one of the other clusters is
analogous.

1. **Download the following files for setting up your cluster account.**

    !!! note

        The following files can be
        [downloaded as an archive here](../_static/cluster_setup_files/cluster_setup_files.tar){:download=cluster_setup_files.tar.}

    With regards to `cluster-setup`, "software scripts" are
    executable files that install and configure software. They are specified
    with the `--software-script` CLI option or `software_scripts` key in the
    `cluster-setup` configuration file.

    - [`install_custom_commands.bash`](../_static/cluster_setup_files/software_scripts/install_custom_commands.bash){:download=install_custom_commands.bash}:
      This script will install a number of useful CLI tools that can be loaded (added to the `PATH` variable)
      via a module. It is accompanied by the [`custom_commands.tar.zst`](../_static/cluster_setup_files/sources/custom_commands.tar.zst){:download=custom_commands.tar.zst}
      archive and [`custom-commands.lua.j2`](../_static/cluster_setup_files/sources/custom_commands.tar.zst){:download=custom-commands.lua.j2}
      template.
    - [`configure_software.bash`](../_static/cluster_setup_files/software_scripts/configure_software.bash){:download=configure_software.bash}:
      This script configures [ASE][ase], [autojob][autojob], and [ccu][ccu] for use by creating suitable
      configuration files. It is accompanied by corresponding templates for
      [ASE](../_static/cluster_setup_files/templates/configuration/ase.ini.j2){:download=ase.ini.j2}
      and [autojob](../_static/cluster_setup_files/templates/configuration/autojob.toml.j2){:download=autojob.toml.j2}
      configuration files and [template scripts](../_static/cluster_setup_files/support_files/autojob_templates.tar.zst){:download=autojob_templates.tar.zst}
      for `autojob`.
    - [`configure_vasp.bash`](../_static/cluster_setup_files/software_scripts/configure_vasp.bash){:download=configure_vasp.bash}: This script configures VASP to be used by ASE
      by copying various support files (pseudopotentials, vDW-DF kernel, and a
      Python script used to call VASP) to a directory.

    !!! note

        Reach out to a group member for instructions of how to get set up with
        the VASP files.

    !!! warning

        If any of above files are downloaded individually (in particular, the
        templates or `custom_commands.tar.zst` archive), be sure to update the
        appropriate path in the corresponding software script.

2. **Copy the files to a DRA cluster.**

    First, collect the files into a directory.

    ```shell
    mkdir cluster_setup_files
    cp FILE_1 FILE_2 FILE_3 ... cluster_setup_files
    ```

    Then, copy the files to a DRA cluster using `scp`.

    ```shell
    scp -r cluster_setup_files <dra_username>@fir.alliancecan.ca
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
        echo "alias $host='ssh -Y $DRA_USER@$host.alliancecan.ca'" >> ~/.zshrc
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
    ls cluster_setup_files
    ```

4. **Create an SSH key and start the `ssh-agent`.**

    In order to interact with GitHub from the cluster, you will need an SSH
    key. You can create an SSH key in the same way that you did on your local
    machine (see the [Basic Setup tutorial](./basic_setup.md)).

    !!! tip

        Don't forget to add the SSH key to your GitHub profile!

    To forego having to type your SSH key every time your require it to
    authenticate to a remote server, start the SSH agent:

    ```shell
    eval $(ssh-agent) && ssh-add path/to/SSH/private/key
    ```

    If you have created an ED25519 SSH key, the path to your private key is
    likely `~/.ssh/ed_25519`.

5. **Create a virtual environment in which to install the `cluster-setup` package.**

    ```shell
    python -m venv cluster_setup_files/.venv && source cluster_setup_files/.venv/bin/activate
    pip install cluster-setup[test]
    ```

    ??? info "Explanation"

        The first command creates a Python virtual environment at `.venv` in the
        current directory and activates it. The final command installs the
        `cluster-setup` package and its `test` extra. Installation
        of the `test` extra enables one to run the test suite prior in
        executing `cluster-setup`.

6. **Run the `cluster-setup` test suite.**

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

7. **Create a configuration file.**

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
      For example, to use the `install_custom_commands.bash`,
      `configure_software.bash` and `configure_vasp.bash` software scripts,
      add something like the following to your configuration file:

      ```toml
      ...
      software_scripts = [
          "path/to/configure_software.bash::::{support_file_home} ~/.config",
          "path/to/configure_vasp.bash::::{support_file_home}",
          "path/to/install_custom_commands.bash:path/to/custom-commands.lua.j2:custom-commands:0.0.1:{software_home}",
      ]
      ...
      ```

8. **Execute `cluster-setup`.**

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

9.  **Verify the setup.**

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

10. **Clean up.**

    The `cluster_setup_files/` directory can now be safely deleted.

    ```shell
    rm -rf cluster_setup_files/
    ```

11. **Rinse and repeat.**

    Log out of the cluster

    ```shell
    exit
    ```

    and repeat steps 2-9 on all clusters on which you would like to run
    calculations.

[cluster-setup]: https://pypi.org/project/cluster-setup/
[ase]: http://ase-lib.org
[autojob]: https://python-autojob.readthedocs.io/en/development/
[ccu]: https://python-comp-chem-utils.readthedocs.io/en/development/
