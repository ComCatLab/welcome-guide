# 🖥 ComCat Lab Welcome Guide (WIP) 🐱

Welcome to the group! This guide introduces several of the key tools and services
that are used in [ComCat Lab][comcat-lab] and examples of their use.

You can view the complete guide in PDF form (comcatlab-welcome-guide.pdf), or,
alternatively, you can [view the guide as a static webpage](#build-webpage).

The best way to learn is by doing! Start by picking a tutorial and trying to work
through it. The tutorials are (relatively) well documented, though, they are not
self-contained. Read over the respective README files to get a cursory overview
of what the tutorial is trying to achieve and then dive in! While you're getting
acquainted with the material and workflows, feel free to bounch back and forth
between the tutorials and the software pages. For more in-depth treatment of the
theory, you are encouraged to consult the "Fundamentals" collection in [ComCat Lab's
Zotero collection][zotero-collection]. There you'll find annotated original references,
more extensive DFT guides, and relevant textbooks. Regarding software, the
documentation pages are generally extensive. Of course, if all else fails, don't
hesitate to reach out to a group member!

## 🌟 Highlights 🌟

This guide has several valuable resources to help you get up to speed:

- **Software Pages**: introductions to the various tools used in ComCat Lab

- **Tutorials**: walkthroughs for several common workflows/softwares

- **Sample Scripts**: example Python and SLURM scripts

- **Troubleshooting**: common issues and their resolutions/work-arounds

## :rocket: Quickstart :rocket:

These instructions assume that you have already completed setup of your SFU
workstation. If you don't have have access to a configured workstation, please out to a
group member to help you set that up, or, if you're feeling up for it, check out
the [local setup repository][local-setup] that will walk you through a comprehensive
setup!

1. Create a virtual environment in the top-level directory of the Welcome Guide folder
(in the same directory as this README).

```bash
python3 -m venv .venv
```

2. Activate your environment.

```bash
source .venv/bin/activate
```

3. Install the required Python dependencies.

```bash
pip install -r tutorial_requirements.txt
```

4. Follow the tutorial!

## ⚒ Build Webpage ⚒

To view the complete guide as a website, you need to install [mkdocs][mkdocs]

```bash
pip install mkdocs
```

then run

```bash
mkdocs build && mkdocs serve
```

[zotero-collection]: https://www.zotero.org/groups/5526800/comcat_lab/library
[local-setup]: https://github.com/ComCatLab/local-setup
[comcat-lab]: https://www.siahrostamilab.com
