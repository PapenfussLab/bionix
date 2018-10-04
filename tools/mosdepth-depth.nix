{ bionix
, nixpkgs
, indexAttrs ? {}
, flags ? null
}:

with nixpkgs;
with lib;

input:

stdenv.mkDerivation {
  name = "mosdepth-depth";
  buildInputs = [ mosdepth ];
  buildCommand = ''
    mkdir $out
    ln -s ${input} input.bam
    ln -s ${bionix.samtools.index indexAttrs input} input.bam.bai
    mosdepth -t $NIX_BUILD_CORES ${optionalString (flags != null) flags} $out/out input.bam
  '';
}
