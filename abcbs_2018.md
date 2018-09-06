We show how Nix, a cross-platform advanced package manager, cleanly solves a number of reproducibility headaches in bioinformatics and computational biology.
Nix can easily create and manage isolated environments, and with our transparent and lightweight extensions can also succinctly describe computational pipelines (workflows), manage their execution in HPC environments or across multiple machines, and produce portable containers (Docker or Singularity images).
We illustrate all our techniques on a typical bioinformatics pipeline.

We compare our approach with the conda software suite. Nix lacks the extensive suite of bioinformatics packages available in Bioconda, but provides a significantly higher degree of reproducibility due to its strong isolation and declarative language.
Moreover, we demonstrate that conda can be used within Nix to leverage Bioconda packages—with some loss of the reproducibility guarantees that a pure Nix solution would entail.


