**Reproducible bioinformatics with Nix**

Justin Bedő, Leon Di Stefano, and Tony Papenfuss

> A challenge for bioinformaticians is to make our computations reproducible — that is, easy to rerun, combine, share, and guaranteed to generate the same results.
> We show how Nix, a next generation cross-platform software deployment system, cleanly overcomes problems usually tackled with a combination of package managers (e.g. conda), containers (Docker, Singularity), and workflow engines (e.g. Toil, Ruffus).
> 
> On its own Nix can be used as a package manager; it can also easily create isolated development environments and export portable containers to share with others.
> We have created a number of transparent and lightweight extensions that enable Nix to succinctly specify bioinformatics analysis environments and pipelines locally, in HPC environments, or in the cloud.
> 
> Nix uses hash-based naming to ensure that what it builds is uniquely-specified, isolation and completeness to ensure that its build processes are deterministic, and a very simple programming language to ensure that the whole system is easy to manage.
> It has an extensive package collection which includes all of CRAN and Bioconductor, and while it lacks Bioconda’s coverage of standalone bioinformatics tools, we show that Bioconda packages can be called from within Nix expressions, with some attendant loss of reproducibility.
> Nix is well supported and general-purpose software that has been in development for over 10 years.
> 
> Here, we will demonstrate how Nix can be used to specify a typical bioinformatics pipeline, and show how it can be executed in whole or in part on an HPC queuing system, or in the cloud.
