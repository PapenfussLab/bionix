# bionix

Bionix is a set of [Nix](http://nixos.org/nix) expressions for specifying and
executing bioinformatics pipelines. It is currently a work in progress, so
documentation is sparse.

## Getting started

Install [Nix](http://nixos.org/nix) and then try `nix-build` in the examples
directory.

## Example pipeline

The examples directory shows a simple pipeline in the call.nix and default.nix
files. You can build it with `nix-build` in the examples directory. The result
will be linked to ./result, which in this example case is an empty VCF file.

## Contact

Please come chat with us at [#bionix:cua0.org](http://matrix.to/#/#bionix:cua0.org).
