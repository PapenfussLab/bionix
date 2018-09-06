We show how Nix, a cross-platform advanced package manager, can cleanly solve a number of reproducibility headaches in bioinformatics and computational biology.
Nix can easily create and manage isolated environments, and with our transparent and extremely lightweight extensions can also describe computational pipelines (workflows), manage their execution in HPC environments, produce containers (Docker or Singularity images) for execution elsewhere, and build in parallel across multiple machines (build farm).
Compared to (Bio)conda, Nix provides a significantly higher degree of reproducibility due to its strong isolation and declarative language as well as remote parallel building capabilities, though lacks the extensive number of bioinformatics packages.
We illustrate our techniques on a typical bioinformatics pipeline.
We show that the entire pipeline including the required software can be specified in a succinct manner and built in parallel either on a local machine directly or via a HPC cluster with a queuing system.
We demonstrate that conda can be used within Nix to leverage the Bioconda repository, with some loss of the reproducibility guarantees a pure Nix solution would entail.
Finally, we discuss how cloud resources can be used to construct a build farm and execute the pipeline.
