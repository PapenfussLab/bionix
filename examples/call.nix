# This is an example pipeline specification to do multi-sample variant calling
# with the Platypus variant caller. Each input is preprocessed by aligning
# against a reference genome (defaults to GRCH38), fixing mate information, and
# marking duplicates. Finally platypus is called over all samples.
{bionix ? import <bionix> {}
,nixpkgs ? import <nixpkgs> {}
,inputs
,ref ? null}:

with bionix;

let
  preprocess = f:
    samtools.markdup {}
      (samtools.sort {}
        (samtools.fixmate {}
          (bwa.align {ref = if ref == null then bionix.ref.grch38.seq else ref;} f)));

  in platypus.call {} (map preprocess inputs)
