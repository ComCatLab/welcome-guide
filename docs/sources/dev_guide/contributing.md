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
- Any page/section marked with "WIP" ("work-in-progress")

## Guidelines for Contributing

When adding new content to tutorials and sample files, please to try
following the existing format. Some general tips:

- Try to limit markdown lines to 80 characters to promote readability. In
  some cases (e.g., URLs, code, etc.), this is not possible.
- Please try to use [fenced code blocks][fenced-code-blocks] to display any
  code that should be copied. This ensures that `mkdocs` formats the code
  correctly and adds copy buttons. This greatly simplifies the user experience.
- Specify a language for fenced code blocks. Often this will simply be `py`,
  `shell`, or `text`. The language is specified inline with the first three
  backticks used to indicate the start of the fenced code block.
- If referring to another document, please consider
  [linking to it directly][linking-to-pages] to make it easier for users to
  navigate.

The following sections have more specific guidelines for adding content to
particular sections.

### Samples

- Add new sample files to an appropriate subdirectory in the top-level
  `samples/` directory and add a corresponding explanation of the sample file
  in the `docs/sources/samples` directory
- Files can be "included" in markdown files using
  ["Snippets Notation"][snippet-notation]. The general syntax is
   `--8<-- ./path/to/file`. `./path/to/file` should be the path to the sample
   file.
- You can add a title to a fenced code block using
  ["Code Block Title Headers"][code-block-headers]. To specify the title,
  add `title=Name of the Title` inline with the code language. Check the
  [Python sample files](../samples/python.md) source code for examples.

### Tutorials

- Specify a "Last Updated" date when writing tutorials. Software is always
  changing, so a tutorial written today may not work tomorrow. Indicating
  when the tutorial was last updated gives users an idea of its currency
  and us an idea of what may need updating.

### Troubleshooting

- Describe the problem as completely as possible. This includes
  listing all relevant software and system information (and versions) and all
  necessary steps to reproduce the problem. If you consulted a resource to
  find the solution, it is also helpful to include links to these resources
  for reference.

### Snippets

- Describe the use-case as clearly as possible
- Annotate your code so that someone who is not so familiar with the language
  will be able to follow along with what each part of the code is trying to
  achieve; if there are any tricky details, this is a good way to highlight
  them!

[git-workflow]: https://www.git-scm.com/book/en/v2/ch00/wfdiag_b
[pull-requests]: https://github.com/ComCatLab/welcome-guide/pulls
[issues]: https://github.com/ComCatLab/welcome-guide/issues
[git-best-practices]: https://about.gitlab.com/topics/version-control/version-control-best-practices/
[fenced-code-blocks]: https://www.markdownguide.org/extended-syntax/#fenced-code-blocks
[linking-to-pages]: https://www.mkdocs.org/user-guide/writing-your-docs/#internal-links
[snippet-notation]: https://facelessuser.github.io/pymdown-extensions/extensions/snippets/#snippets-notation
[code-block-headers]: https://facelessuser.github.io/pymdown-extensions/extensions/superfences/#code-block-title-headers
