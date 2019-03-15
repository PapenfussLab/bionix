{ bionix
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, alignerStreaming ? false
, flags ? null
, config ? null
, heapSize ? "2G"

}:

with bionix;
with lib;
with types;

input:

let
  ref = matchFiletype "gridss-softClipsToSplitReads" { bam = x: x.ref; } input;
in

stage rec {
  name = "gridss-softClipsToSplitReads";
  buildInputs = with pkgs; [ jre bwa ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    java -Xmx${heapSize} -Dsamjdk.create_index=false \
      -cp ${bionix.gridss.jar} gridss.SoftClipsToSplitReads \
			REFERENCE_SEQUENCE=ref.fa \
			I=${input} \
			O=$out \
      ${optionalString alignerStreaming "ALIGNER_STREAMING=true"} \
      ${optionalString (config != null) ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)} \
			WORKER_THREADS=$NIX_BUILD_CORES
    '';
  passthru.filetype = filetype.bam { ref = ref; sorting = sort.none {}; };
  passthru.multicore = true;
}
