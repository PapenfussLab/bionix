{ bionix
, nixpkgs
, flags ? null
}:

input:

with nixpkgs;
with lib;

stdenv.mkDerivation {
  name = "samtools-index";
  buildInputs = [ samtools ];
  buildCommand = ''
    ln -s ${input} input.bam
    samtools index -@ $NIX_BUILD_CORES ${optionalString (flags != null) flags} input.bam
    cp input.bam.bai $out
  '';
}
