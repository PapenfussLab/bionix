# This is an example tumour-normal calling pipeline using strelka
{ bionix ? import ./../.. {}
, normal
, tumour
, ref
}:

with bionix;
with lib;

let
  input = mapAttrs (_: fetchFastQGZ);

  preprocess = pipe [
    input
    (bwa.align { ref = fetchFastA ref; })
    (samtools.fixmate {})
    (samtools.sort {})
    (samtools.markdup {})
  ];

in linkOutputs {
  strelka = strelka.callSomatic {} {normal = preprocess normal; tumour = preprocess tumour;};
  "normal.bam" = preprocess normal;
  "tumour.bam" = preprocess tumour;
}
