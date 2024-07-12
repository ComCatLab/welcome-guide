# How to Contribute

This page highlights sections of the guide that need work
and explains how to report issues/suggest new features.

## Development Workflow

Assuming that you have [set up your environment](setup):

1. Create and checkout a branch for your local changes.

    ```bash
    git checkout -b name-of-feature-or-change
    ```

    !!! reminder

        It is really important that you do not push your changes directly to the `main`
        branch. In the case that someone makes a mistake with `git`, it is much easier
        to correct the issue on a single fork (these can be deleted and re-forked
        relatively painlessly) than for the whole project. See
        [Integration-Manager Workflow][git-workflow] for a more detailed explanation
        of the workflow.)

2. Make your changes (e.g., add/change/remove files).

3. Commit your changes.

    ```bash
    git commit -S -m "I made a change"
    ```

4. Push your changes to **your** remote.

    ```bash
    git push origin name-of-feature-or-change
    ```

5. [Create a pull request][pull-requests] for your change with an appropriate template.

## Where to Contribute

- Writing GitHub templates for issues/MRs
- Configuring GH actions
- Writing software pages
- Creating tutorials for new calculations
- Uploading sample files
- Explaining solutions to technical issues
- Compiling useful links

[git-workflow]: https://www.git-scm.com/book/en/v2/ch00/wfdiag_b
[pull-requests]: https://github.com/ComCatLab/welcome-guide/pulls
