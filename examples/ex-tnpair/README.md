This example is a tumour-normal processing workflow applied to a
publically available whole-genome sequencing (WGS) dataset. As this
example uses WGS data, a large amount of data will be downloaded.
Furthermore, a large amount of space will be required to build the final
products.

# Building on local machine

Run `nix build -I bionix=../..` in this directory.

# Building via HPC (slurm)

Run `nix build -f cluster.nix` to build on slurm. Note that Nix must be
configured such that the temporary build directories are created on
shared storage. You may need to adjust the resource limits specified in
cluster.nix to suit your particular cluster hardware.
