{ bionix
, nixpkgs
, dictIndexAttrs ? {}
, faidxAttrs ? {}
, flags ? null
, unmappedReads ? false
, minClipLength ? 5
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
    mkdir $out
    java -Dsamjdk.create_index=true \
      -cp ${bionix.gridss.jar} gridss.ExtractSVReads \
      REFERENCE_SEQUENCE=ref.fa \
      I=input.bam \
      O=$out/input.sv.bam \
      METRICS_OUTPUT=$out/input.sv_metrics \
      INSERT_SIZE_METRICS=$out/input.insert_size_metrics \
      UNMAPPED_READS=${if unmappedReads then "true" else "false"} \
      MIN_CLIP_LENGTH=${toString minClipLength}
  '';
}
