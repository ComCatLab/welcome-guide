---
hide:
  - toc
---

# 🖥 ComCat Lab Welcome Guide 🐱

Welcome to the group! This guide introduces several of the key tools and services
that are used in [ComCat Lab][comcat-lab] and examples of their use.

You can view the complete guide in PDF form (`comcatlab-welcome-guide.pdf`), or,
alternatively, you can [view the guide as a static webpage](dev_guide/docs.md#building-the-webpage).

The best way to learn is by doing! Start by picking a tutorial and trying to work
through it. The tutorials are (relatively) well documented, though, they are not
self-contained. While you're getting acquainted with the material and workflows,
feel free to bounce back and forth between the tutorials and the software pages.
For more in-depth treatment of the theory, you are encouraged to consult the
"Fundamentals" collection in [ComCat Lab's Zotero collection][zotero-collection].
There you'll find annotated original references, more extensive DFT guides, and references
to relevant textbooks. Regarding software, the documentation pages are generally extensive.
Of course, if all else fails, don't hesitate to reach out to a group member!

## 🌟 Highlights 🌟

This guide has several valuable resources to help you get up to speed:

- [**On-Boarding**](onboarding/index.md): your first steps to getting started in ComCat Lab

- [**Workflows**](nutshell.md): an overview of frequent workflows used with links to relevant software pages

- [**Software Pages**](software_pages.md): introductions to the various tools used in ComCat Lab

- [**Tutorials**](tutorials/index.md): walkthroughs for several common workflows/softwares

- **Sample Scripts**: example [Python](samples/python.md) and
  [SLURM](samples/slurm.md) scripts

- [**Tips**](resources/troubleshooting.md): common issues and their resolutions/work-arounds

## Prerequisites

- [**Set up your local machine**][local-setup]: this guide assumes that your local machine
  is setup up with the necessary software already installed (e.g., Hatch, Python, etc.);
  follow the steps outlined in [local-setup][local-setup] to satisfy this prerequisite

- [**Set up your cluster account**][cluster-setup]: some of the tutorials in the guide
  require that you connect to the remote clusters provided by the Alliance. For these
  tutorials, you will need a [valid CCDB account](tutorials/ccdb.md) and to
  [set up your remote cluster environment][cluster-setup]

[comcat-lab]: https://www.siahrostamilab.com
[zotero-collection]: https://www.zotero.org/groups/5526800/comcat_lab/library
[local-setup]: https://github.com/ComCatLab/local-setup
[cluster-setup]: https://github.com/ComCatLab/cluster-setup
