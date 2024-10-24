# Building the Documentation

This page describes how the documentation pages of the Welcome Guide
are built.

## Build Tools

This guide uses `mkdocs` to convert the Markdown files to HTML files
and build the webpage. All the configuration for using `mkdocs`
is contained within the `mkdocs.yml` configuration file in the project root.
For a detailed description, please see the [mkdocs documentation][mkdocs].
For theme options, please see the documentation for
[material for mkdocs][material].

## Building the Webpage

To build the webpage,

1. Activate the development environment.

    If you are using [Hatch][install-hatch]:

    ```bash
    hatch shell
    ```

    If you created a virtual environment with `venv` (e.g.,
    `python -m venv .venv`):

    ```bash
    source .venv/bin/activate
    ```

    Note that at the very least, you should have the `docs` extra installed
    (e.g., `pip install '.[docs]'`).

2. Build the HTML files.

    ```bash
    mkdocs build
    ```

    The output should be something like this:

    ```bash
    INFO    -  Cleaning site directory
    INFO    -  Building documentation to directory: /Users/ugo/Projects/nwt/welcome-guide/docs/site
    INFO    -  The following pages exist in the docs directory, but are not included in the "nav" configuration:
                 - resources/conferences.md
                 - resources/links.md
                 - resources/references.md
                 - resources/social.md
                 - software_pages/index.md
    INFO    -  Documentation built in 0.19 seconds
    ```

3. Serve the webpage.

    ```bash
    mkdocs serve
    ```

    The output should be something like this:

    ```bash
    INFO    -  Building documentation...
    INFO    -  Cleaning site directory
    INFO    -  The following pages exist in the docs directory, but are not included in the "nav" configuration:
                 - resources/conferences.md
                 - resources/links.md
                 - resources/references.md
                 - resources/social.md
                 - software_pages/index.md
    INFO    -  Documentation built in 0.16 seconds
    INFO    -  [16:49:00] Watching paths for changes: 'docs/sources', 'mkdocs.yml'
    INFO    -  [16:49:00] Serving on http://127.0.0.1:8000/
    ```

4. Enter the IP address into your browser (here, `127.0.0.1:8000/`)

## Building the PDF Version of the Guide

To build the PDF version of the guide, you will need to install
[Pandoc][pandoc] and a [LaTeX][latex] engine. You will then need to remove all
emoji symbols from the markdown source documents. Finally, from the project
root, run:

[//]: # (Update this command whenever new files are added to the docs/source/ directory)

```bash
cd docs/sources || exit; pandoc --file-scope -s -o ../../comcat-lab-welcome-guide.pdf -f markdown -t pdf index.md quickstart.md nutshell.md software_pages.md samples/{index,python,bash,slurm}.md tutorials/{index,*}.md resources/*.md dev_guide/*.md; cd ../../ || exit
```

If you have hatch configured, you can run:

```bash
hatch run docs:build-pdf
```

This command will collate the markdown source files into a single PDF **in the
order that their filenames** are passed to the `pandoc` command. The PDF will
be saved to the project root in a file called `comcat-lab-welcome-guide.pdf`.

!!! warning

    As of 2024/10/23, the formatting for this is still problematic. Code
    formatting is poor, included files (e.g., in `samples/`) are not rendered,
    and symbols result in errors.

[mkdocs]: https://www.mkdocs.org/user-guide/
[material]: https://squidfunk.github.io/mkdocs-material/setup/
[pandoc]: https://pandoc.org/installing.html
[latex]: https://www.latex-project.org/get/
[install-hatch]: https://hatch.pypa.io/latest/install/
