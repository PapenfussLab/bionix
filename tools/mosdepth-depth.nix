{ stdenv
, lib
, callPackage
, mosdepth
, index ? callPackage ./samtools-index.nix {}
, flags ? null}:

with lib;

input:

stdenv.mkDerivation {
  name = "mosdepth-depth";
  buildInputs = [ mosdepth ];
  buildCommand = ''
    mkdir $out
    ln -s ${input} input.bam
    ln -s ${index input} input.bam.bai
    mosdepth -t $NIX_BUILD_CORES ${optionalString (flags != null) flags} $out/out input.bam
  '';
}
