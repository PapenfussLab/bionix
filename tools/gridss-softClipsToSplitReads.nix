{ bionix
, nixpkgs
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, alignerStreaming ? false
, flags ? null
}:

with nixpkgs;
with lib;
with bionix.types;

input:

let
  ref = matchFiletype "gridss-softClipsToSplitReads" { bam = x: x.ref; } input;
in

assert (matchFileSorting "gridss-softClipsToSplitReads" { name = _: true; } input);

stdenv.mkDerivation rec {
  name = "gridss-softClipsToSplitReads";
  buildInputs = [ jre ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    java -Xmx2G -Dsamjdk.create_index=false \
      -cp ${bionix.gridss.jar} gridss.SoftClipsToSplitReads \
			REFERENCE_SEQUENCE=ref.fa \
			I=${input} \
			O=$out \
      ${optionalString alignerStreaming "ALIGNER_STREAMING=true"} \
			WORKER_THREADS=$NIX_BUILD_CORES
    '';
  passthru.filetype =
    if alignerStreaming then
      filetype.bam { ref = ref; sort = sorting.none {};  }
    else
      input.filetype;
}