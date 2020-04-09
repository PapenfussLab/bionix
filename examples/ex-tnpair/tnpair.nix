{bionix ? import <bionix> {}, pair, fetch}:

with bionix;
with lib;
with types;

with minimap2;
with samtools;
with snpeff;

let
  preprocess = s: pipe s [
    fetch
    (align { preset = "sr"; ref = ref.grch38.seq; flags = "-R'@RG\\tID:${s.type}\\tSM:${s.type}'"; })
    (fixmate {})
    (sort { })
    (markdup { })
  ];

  dropErrors = input: stage {
    name = "drop-errors";
    buildCommand = ''
      grep -v "ERROR_" ${input} > $out
    '';
    passthru.filetype = input.filetype;
  };

  bams = mapAttrs (_: preprocess) pair;

  variants = let
    somatic = strelka.callSomatic { } bams; in mapAttrs (_: flip pipe [
      (compression.uncompress { })
      (snpeff.annotate { db = ref.grch38.snpeff.db; })
      dropErrors
      (snpeff.dbnsfp { dbnsfp = ref.grch38.snpeff.dbnsfp; })
    ]) {
      "snvs.vcf" = somatic.snvs;
      "indels.vcf" = somatic.snvs;
      "germline.vcf" = strelka.call { } [bams.normal];
    };

  cnvs = cnvkit.callCNV { } { normals = [ bams.normal ]; tumours = [ bams.tumour ]; };

in linkOutputs {
  inherit variants;
  alignments = linkOutputs (mapAttrs' (n: nameValuePair (n + ".bam")) bams);
  cnvkit = cnvs;
}
