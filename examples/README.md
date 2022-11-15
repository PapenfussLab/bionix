# Bionix examples

This directory has a few example workflows in bionix along with example
data. A basic workflow is defined in `call.nix`, and an example of
applying it to the sample data is in `default.nix`. To build the
`default.nix` workflow, run ```nix build -I bionix=../ -f .``` from this directory.

## NextFlow and WDL translations

The directories `ex-nextflow` and `ex-wdl` contain translated examples
from the NextFlow and WDL documentation respectively.

The NextFlow translated example does not come with example data. It can be built with
```
nix build -f nextflow-example1.nix -I bionix=../.. --arg input /path/to/sample.fa
```

The WDL example requires no extra data and can be built with
```
nix build -f wdl-scatter-gather.nix -I bionix=../..
```

## Example script wrapper

`ex-tnpair` contains a shell script based example on how a front-end for
users might be constructed. It is a simple tumour-normal somatic calling
workflow using the Strelka variant caller. The script accepts a
reference fasta along with paired normal and tumor fastq files and
performs alignment, preprocessing, and variant calling with
[`strelka`](https://github.com/Illumina/strelka).
