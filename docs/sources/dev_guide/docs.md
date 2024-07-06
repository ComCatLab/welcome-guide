# Building the Documentation

This page describes how the documentation pages of the Welcome Guide
are built.

## Build Tools

This guide uses `mkdocs` to convert the Markdown files to HTML files
and build the webpage. All the configuration for using mkdocs
is contained within the `mkdocs.yml` configuration file in the project root.
For a detailed description, please see the [mkdocs documentation][mkdocs].

## Building the Webpage

To build the webpage, ensure that `mkdocs` and the necessary extensions specified in
the `mkdocs.yml` configuration file are installed in the current Python environment.
For convenience, the `docs` Hatch environment is maintained to contain all necessary
Python dependencies for building the documentation. To build to webpage, run

```bash
>>> mkdocs build
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

from within a suitable Python environment.

## Serving the Webpage

To view a live version of the Welcome Guide as a webpage, ensure that `mkdocs` and the
necessary extensions specified in the `mkdocs.yml` configuration file are installed in
the current Python environment. For convenience, the `docs` Hatch environment is maintained
to contain all necessary Python dependencies for serving the documentation. To serve the
webpage, run

```bash
>>> mkdocs serve
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

[mkdocs]: https://www.mkdocs.org/user-guide/
