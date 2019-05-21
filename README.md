<h1 align="center"> BioNix </h1>

BioNix is a tool for reproducible bioinformatics that unifies workflow engines, package managers, and containers.
It is implemented as a lightweight library on top of the [Nix](https://nixos.org/nix/) deployment system.

BioNix is currently a work in progress, so documentation is sparse.

## Getting started

Install [Nix](http://nixos.org/nix):

```{sh}
curl https://nixos.org/nix/install | sh
```
To run a sample pipeline, clone this project and run `nix-build` in the `/examples` directory:

```{sh}
$ git clone https://github.com/PapenfussLab/bionix
$ cd examples
$ nix-build
```

The sample pipeline performs variant calling using [`platypus`](https://github.com/andyrimmer/Platypus), alignment using [`bwa mem`](https://github.com/lh3/bwa), and preprocessing using [`samtools`](http://www.htslib.org/). 
BioNix will download or build all of the necessary software and create a soft link (`result`) to the workflow output.

Next, check out the code:

-   The pipeline itself is specified in `examples/call.nix` and `examples/default.nix`.
-   The BioNix wrapper for `platypus` is in `tools/platypus-callVariants.nix`.
-   The software package for `platypus` can be found in [nixpkgs](https://github.com/NixOS/nixpkgs/blob/master/pkgs/applications/science/biology/platypus/default.nix).

BioNix pipelines can be easily wrapped in shell scripts: see `examples/ex-tnpair/tnpair` for an example script that accepts a reference fasta, along with paired normal and tumor fastq files, and performs alignment, preprocessing, and variant calling with [`strelka`](https://github.com/Illumina/strelka).

Writing your own pipelines requires some familiarity with the Nix programming language and deployment system. Good introductions can be found [here](https://learnxinyminutes.com/docs/nix/) and [here](https://ebzzry.io/en/nix/).

We have successfully run BioNix pipelines in a zero-install manner (using a [statically linked binary](https://matthewbauer.us/blog/static-nix.html) and [user namespaces](https://www.redhat.com/en/blog/whats-next-containers-user-namespaces)), but this feature is currently unstable. Stay tuned!

## Contact

Please come chat with us at [#bionix:cua0.org](http://matrix.to/#/#bionix:cua0.org).
