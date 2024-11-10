# Snippets

This page contains useful snippets for performing common tasks.

## Modulefiles

### Create a modulefile from a script

[Lmod](../software_pages.md#lmod) provides a handy utility
[`sh_to_modulefile`][sh_to_modulefile] for creating modulefiles from scripts.
`sh_to_modulefile` records the initial environment, runs the script, and
compares the final environment to determine the changes. It then converts these
changes to Lua commands and prints them to the terminal. Given a shell script,
`my_script.sh`, one can run the following command to generate an Lmod-compatible
modulefile written in Lua.

```shell
$LMOD_DIR/sh_to_modulefile my_script.sh > my_script.lua
```

The command `$LMOD_DIR/sh_to_modulefile` calls the utility by using the
Lmod-defined environment variable `LMOD_DIR`, which points to the directory
in which Lmod is installed. `my_script.sh` is the name of the script to be
converted into a modulefile. The output of the command is redirected (`>`)
to the file `my_script.lua`.

[sh_to_modulefile]: https://lmod.readthedocs.io/en/latest/260_sh_to_modulefile.html

## SLURM

### Attaching to Existing Jobs

Say, you have a job running on Cedar (this could be an interactive job launched
with `salloc` or a job submitted with `sbatch`) and you would like to run some
commands on the compute node with the job's resources (cores, memory, GPU, etc.).
If the job ID for the job is JOBID, you can "attach" to the job with the command:

```shell
srun --pty --overlap --jobid JOBID bash
```

This snippet calls the [`srun` command of SLURM][srun]. The `--pty` option
launches a task in "pseudo terminal mode", allowing for you to specify commands
from the command line as if you were running a terminal session on the compute
node. The `--overlap` option allows your terminal session to share resources
with the running job. The `--jobid` option specifies the SLURM job ID to
attach to.

!!! warning

    Since the commands that you run in this terminal will share the resources
    of the running job, if you exceed the memory request for the job, you can
    cause your original job to fail due to an out-of-memory (OOM) error.

[srun]: https://slurm.schedmd.com/srun.html
