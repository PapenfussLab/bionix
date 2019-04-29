# This is an example tumour-normal calling pipeline using strelka
{ bionix ? import <bionix> {}
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

in linkDrv [
  (ln (strelka.call {} {normal = preprocess normal; tumour = preprocess tumour;}) "strelka")
  (ln (preprocess normal) "normal.bam")
  (ln (preprocess tumour) "tumour.bam")
]
