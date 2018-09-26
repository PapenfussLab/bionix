{ stdenv
, callPackage
, lib
, samtools
}:

input:

with lib;

stdenv.mkDerivation {
  name = "samtools-index";
  buildInputs = [ samtools ];
  buildCommand = "samtools flagstat -@ $NIX_BUILD_CORES ${input} > $out";
}
