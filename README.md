# bionix

Bionix is a set of [Nix](http://nixos.org/nix) expressions for specifying and
executing bioinformatics pipelines. It is currently a work in progress, so
documentation is sparse.

## Getting started

Install [Nix](http://nixos.org/nix) and then try `nix-build` in the examples
directory.

## Example pipelines

The examples directory shows a simple pipeline in the call.nix and default.nix
files. You can build it with `nix-build` in the examples directory. The result
will be linked to ./result, which in this example case is an empty VCF file.

There is another tumour-normal pipeline calling example consisting of `tnpair`
and `tnpair.nix`. In this case, the shell script tnpair accepts a reference and
two fastq files to run the pipeline defined in `tnpair.nix` on.

## Contact

Please come chat with us at [#bionix:cua0.org](http://matrix.to/#/#bionix:cua0.org).
