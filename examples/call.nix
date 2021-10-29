# This is an example pipeline specification to do multi-sample variant calling
# with the Platypus variant caller. Each input is preprocessed by aligning
# against a reference genome (defaults to GRCH38), fixing mate information, and
# marking duplicates. Finally platypus is called over all samples.
{ bionix ? import <bionix> { }
, inputs
, ref ? bionix.ref.grch38.seq
}:

with bionix;
with lib;

let
  preprocess = flip pipe [
    (bwa.align { inherit ref; })
    (samtools.sort { nameSort = true; })
    (samtools.fixmate { })
    (samtools.sort { })
    (samtools.markdup { })
  ];

in
platypus.call { } (map preprocess inputs)
