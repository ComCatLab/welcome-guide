# Contributing to the Welcome Guide (WIP)

First off, thanks for taking the time to contribute! ❤️

All types of contributions are encouraged and valued. See the
[Table of Contents](#table-of-contents) for different ways to help
and details about how this project handles them. Please make sure to read
the relevant section before making your contribution. It will make it a lot
easier for us maintainers and smooth out the experience for all involved. The community
looks forward to your contributions. 🎉

## Table of Contents

- [Contributing to the Welcome Guide (WIP)](#contributing-to-the-welcome-guide-wip)
  - [Table of Contents](#table-of-contents)
  - [Code of Conduct (WIP)](#code-of-conduct-wip)
  - [Getting Started](#getting-started)
  - [Reporting Bugs](#reporting-bugs)
  - [Got a Question?](#got-a-question)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Where to Contribute](#where-to-contribute)

## Code of Conduct (WIP)

- see quacc

## Getting Started

1. Contact an admin to add your SSH key to get access

2. Clone the repository

    ```bash
        git clone git@github.com:ComCatLab/welcome-guide.git
    ```

3. Download [Hatch][hatch]

    ```bash
    brew install hatch
    ```

4. Create the default environment and install the project with Hatch

    ```bash
    hatch shell
    ```

5. Install pre-commit hooks

    pre-commit install --hook-type pre-push --hook-type pre-commit

## Reporting Bugs

Before reporting a bug, please perform the following checks:

1. Search the existing issues to see if someone has already raised the same issue.

2. Ensure that your problem persists on the latest version of the code.

3. Create an issue using an issue template.

## Got a Question?

- open an issue

## Suggesting Enhancements

- open a PR
- hotfixes target main
- new features target development

## Where to Contribute

- add new resources to the collection (and please annotate!)

[hatch]: https://hatch.pypa.io/latest
