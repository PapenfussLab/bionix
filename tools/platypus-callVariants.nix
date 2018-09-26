{ stdenv
, callPackage
, lib
, platypus
, ref
, index ? callPackage ./samtools-faidx.nix {}
, bamIndex ? callPackage ./samtools-index.nix {}
, flags ? null
}:

inputs:

with lib;

let filename = path: last (splitString "/" path);

in stdenv.mkDerivation {
  name = "platypus";
  buildInputs = [ platypus ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${index ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bamIndex p} ${filename p}.bai") inputs}
    ls -l
    platypus callVariants \
      --nCPU=$NIX_BUILD_CORES \
      --refFile=ref.fa \
      ${optionalString (flags != null) flags} \
      -o $out \
      --bamFiles=${concatMapStringsSep "," (p: "${filename p}.bam") inputs}
  '';
}
