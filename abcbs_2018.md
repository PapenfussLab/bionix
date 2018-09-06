**Reproducible bioinformatics with Nix**

Justin Bedő, Leon Di Stefano, and Tony Papenfuss

> A cornerstone of science is reproducibility, the ability to independently verify experimental research. For bioinformatics to support scientific reproducibility, the computational portion of a research project has to be well specified and recomputable. However, it is difficult to guarantee reproducibility for a bioinformatics pipeline, in part due to the large number of software invoked, their complicated interactions, and the size of our data. Recent approaches such as containerisation does not solve this problem as it simply shifts the difficulty to managing containers instead of managing software artifacts. Furthermore, the execution of a pipeline is usually disjoint from the container construction, adding further management difficulties.

> We show how Nix, a next generation cross-platform software deployment system, can cleanly solve a number of reproducibility headaches in bioinformatics and computational biology.
> Nix uses hash-based naming to ensure that its builds are uniquely specified, isolation and completeness to ensure that they are deterministic, and a simple programming language to ensure that they are easily manageable.
> Nix is well supported and mature software with a large community that has been in development for over 10 years.
> 
> With our transparent and lightweight extensions Nix succinctly describe computational pipelines, manage their execution in HPC environments or in parallel across a collection of machines, and produce containers (e.g., Docker, Singularity) to share with others.
> Nix has an extensive package collection that includes the whole of CRAN and Bioconductor, which can be leveraged in our pipelines.
> While Nix lacks Bioconda's coverage of standalone bioinformatics tools, we show that Bioconda can be used within Nix expressions, with some attendant loss of reproducibility.
>
> In our talk we will use Nix to specify a typical bioinformatics pipeline, and show how it can be executed in whole or in part on an HPC queuing system.
