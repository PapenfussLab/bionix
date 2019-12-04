{ bionix
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

inputs:

with bionix;
with lib;
with types;


let
  filename = path: last (splitString "/" path);
  getref = f: matchFiletype "platypus-callVariants" { bam = {ref, ...}: ref; } f;
  refs = map getref inputs;
  ref = head refs;
in

assert (length (unique refs) == 1);

stage {
  name = "platypus";
  buildInputs = with pkgs; [ platypus ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx indexAttrs ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bionix.samtools.index bamIndexAttrs p} ${filename p}.bai") inputs}
    platypus callVariants \
      --nCPU=$NIX_BUILD_CORES \
      --refFile=ref.fa \
      ${optionalString (flags != null) flags} \
      -o $out \
      --bamFiles=${concatMapStringsSep "," (p: "${filename p}.bam") inputs}

    # Remove timestamps from output
    sed -i '/^##fileDate/d' $out
  '';
  passthru.filetype = filetype.vcf {ref = ref;};
  passthru.multicore = true;
}
