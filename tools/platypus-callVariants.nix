{ bionix
, nixpkgs
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

inputs:

with nixpkgs;
with lib;
with bionix.types;


let
  filename = path: last (splitString "/" path);
  getref = f: matchFiletype "platypus-callVariants" { bam = r: r; } f;
  refs = map getref inputs;
  ref = head refs;
in

assert (length (unique refs) == 1);

stdenv.mkDerivation {
  name = "platypus";
  buildInputs = [ platypus ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx indexAttrs ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bionix.samtools.index bamIndexAttrs p} ${filename p}.bai") inputs}
    ls -l
    platypus callVariants \
      --nCPU=$NIX_BUILD_CORES \
      --refFile=ref.fa \
      ${optionalString (flags != null) flags} \
      -o $out \
      --bamFiles=${concatMapStringsSep "," (p: "${filename p}.bam") inputs}
  '';
  passthru.filetype = filetype.vcf {ref = ref;};
}
