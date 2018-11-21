{ bionix
, nixpkgs
, blacklist ? null
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, flags ? null
, config ? null
}:

with nixpkgs;
with lib;
with bionix.types;

input:

let
  ref = matchFiletype "gridss-computeSamTags" { bam = x: x.ref; } input;
  sorted = matchFileSorting "gridss-computeSamTags" { name = _: true; } input;
in

assert(sorted);

stdenv.mkDerivation rec {
  name = "gridss-computeSamTags";
  buildInputs = [ jre ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    java -Xmx1G \
      -Dsamjdk.create_index=false \
			-cp ${bionix.gridss.jar} gridss.ComputeSamTags \
			REFERENCE_SEQUENCE=ref.fa \
			WORKING_DIR=$TMP_DIR \
			TMP_DIR=$TMP_DIR \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + gridssConfig config)} \
			I=${input} \
			O=$out \
      AS=true
  '';
  passthru.filetype = input.filetype;
}
