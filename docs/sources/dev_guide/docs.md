# Building the Documentation

This page describes how the documentation pages of the Welcome Guide
are built.

## Build Tools

This guide uses `mkdocs` to convert the Markdown files to HTML files
and build the webpage. All the configuration for using `mkdocs`
is contained within the `mkdocs.yml` configuration file in the project root.
For a detailed description, please see the [mkdocs documentation][mkdocs].
For theme options, please see the documentation for [material for mkdocs][material].

## Building the Webpage

To build the webpage, ensure that `mkdocs` and the necessary extensions specified in
the `mkdocs.yml` configuration file are installed in the activate Python environment.
For convenience, the `docs` Hatch environment is maintained to contain all necessary
Python dependencies for building the documentation. To build the webpage, run

```bash
mkdocs build
```

from within a suitable Python environment. The output should be something like this:

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

## Serving the Webpage

To view a live version of the Welcome Guide as a webpage, ensure that `mkdocs` and the
necessary extensions specified in the `mkdocs.yml` configuration file are installed in
the activate Python environment. For convenience, the `docs` Hatch environment is maintained
to contain all necessary Python dependencies for serving the documentation. To serve the
webpage, run

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

Now, enter the IP address shown, (here, `127.0.0.1:8000/`), into your browser to view the
webpage.

## Building the PDF Version of the Guide

To build the PDF version of the guide, you will need to install [Pandoc][pandoc] and a [LaTeX][latex]
engine. Then run

```bash
cd docs/sources && pandoc --file-scope -s -o ../../comcat-lab-welcome-guide.pdf -f markdown -t pdf index.md quickstart.md workflows.md software_pages.md samples/*.md tutorials/*.md resources/{troubleshooting.md,links.md} dev_guide/* && cd ../../
```

[mkdocs]: https://www.mkdocs.org/user-guide/
[material]: https://squidfunk.github.io/mkdocs-material/setup/
[pandoc]: https://pandoc.org/installing.html
[latex]: https://www.latex-project.org/get/
