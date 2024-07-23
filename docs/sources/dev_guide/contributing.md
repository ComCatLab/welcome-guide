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

3. Commit your changes. Each commit should include a brief (50 characters or less)
   description of the changes. If you can't describe your changes in 50 characters
   or less, you should probably break the commit into smaller pieces!
   [This][git-best-practices] is a great guide to good commits.

    ```bash
    git commit -S -m "Added VASP relaxation python script"
    ```

4. Push your changes to **your** remote.

    ```bash
    git push origin name-of-feature-or-change
    ```

5. [Create a pull request][pull-requests] with a summary of your changes using
   the pull request template. Please remember to complete the checklist!

## Issues

Maybe you notice an error in the documentation, or maybe, a new version
of a particular software causes one of the tutorials to no longer work.
This is a perfect chance to [raise an issue][issues]! Please use an issue template
to describe the problem and try to include the following information, where
relevant:

- **System**: operating system, software version (e.g., ASE 3.23.0 or VASP 5.4.4)
- **Intended behaviour**: what should be happening?
- **Problem**: what's actually happening?
- **Minimum Reproducible Rxample (MRE)**: once you have an example of the problem,
  try to remove as much as possible from your example while still retaining the error

## Features

If you think that adding something would improve the guide, suggest it by
[submitting an issue][issues]. Please explain what the change would be and how it
would benefit the guide. From there, you can [create a pull request][pull-requests]
with the change!

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
[issues]: https://github.com/ComCatLab/welcome-guide/issues
[git-best-practices]: https://about.gitlab.com/topics/version-control/version-control-best-practices/
