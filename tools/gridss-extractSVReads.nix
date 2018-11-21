{ bionix
, nixpkgs
, dictIndexAttrs ? {}
, faidxAttrs ? {}
, flags ? null
, unmappedReads ? false
, minClipLength ? 5
, config ? null
}:

with nixpkgs;
with lib;
with bionix.types;

input:

let
  ref = matchFiletype "gridss-extractSVReads" { bam = x: x.ref; } input;
in


stdenv.mkDerivation rec {
  name = "gridss-extractSVReads";
  buildInputs = [ jre R ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    ln -s ${bionix.samtools.dict dictIndexAttrs ref} ref.fa.dict
    ln -s ${input} input.bam
    for f in ${bionix.gridss.collectMetrics collectMetricsAttrs input}/* ; do
      ln -s $f
    done
    java -Dsamjdk.create_index=true \
      -cp ${bionix.gridss.jar} gridss.ExtractSVReads \
      REFERENCE_SEQUENCE=ref.fa \
      I=input.bam \
      O=$out \
      UNMAPPED_READS=${if unmappedReads then "true" else "false"} \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + gridssConfig config)} \
      MIN_CLIP_LENGTH=${toString minClipLength}
  '';
  passthru.filetype = input.filetype;
}
