{ bionix
, nixpkgs
}:

input:

with nixpkgs;
with lib;

stdenv.mkDerivation {
  name = "samtools-index";
  buildInputs = [ samtools ];
  buildCommand = "samtools flagstat -@ $NIX_BUILD_CORES ${input} > $out";
}
