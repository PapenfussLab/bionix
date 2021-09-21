{ bionix, bwaIndexAttrs ? { }, faidxAttrs ? { }, flags ? null, config ? null
, heapSize ? "1G" }:

with bionix;
with lib;
with types;

input:

let
  ref = matchFiletype "gridss-computeSamTags" { bam = x: x.ref; } input;
  sorted = matchFileSorting "gridss-computeSamTags" { name = _: true; } input;

in assert (sorted);

stage rec {
  name = "gridss-computeSamTags";
  buildInputs = with pkgs; [ jre ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    java -Xmx${heapSize} \
      -Dsamjdk.create_index=false \
      -cp ${bionix.gridss.jar} gridss.ComputeSamTags \
      VERBOSITY=WARNING \
      WORKER_THREADS=$NIX_BUILD_CORES \
      REFERENCE_SEQUENCE=ref.fa \
      WORKING_DIR=$TMP_DIR \
      TMP_DIR=$TMP_DIR \
      ${
        optionalString (config != null)
        ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)
      } \
      I=${input} \
      O=$out \
      AS=true
  '';
  passthru.filetype = input.filetype;
  passthru.multicore = true;
}
