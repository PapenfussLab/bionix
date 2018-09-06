**Reproducible bioinformatics with Nix**

Justin Bedő, Leon Di Stefano, and Tony Papenfuss

> We show how Nix, a next generation cross-platform software deployment system, cleanly solves a number of reproducibility headaches in bioinformatics and computational biology.
> First developed in 2006, Nix uses hash-based naming to ensure that its builds are uniquely-specified, isolation and completeness to ensure that they are deterministic, and a very simple programming language to ensure that they are easily-manageable.
> 
> Nix, like tools such as (mini)conda, can straightforwardly create and manage isolated environments, and with our transparent and lightweight extensions it can also succinctly describe computational pipelines, manage their execution in HPC environments or in parallel across a collection of machines, and produce portable containers (Docker or Singularity images) to share with others.
> Nix has an extensive package collection which includes the whole of CRAN and Bioconductor, and while it lacks Bioconda's coverage of standalone bioinformatics tools, we show that Bioconda can be used within Nix expressions (with some attendant loss of reproducibility).
>
> In our talk we will use Nix to specify a typical bioinformatics pipeline, and show how it can be executed in whole or in part on an HPC queuing system.
