{ bionix
, nixpkgs
, flags ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-index" { bam = _: true; } input);
assert (matchFileSorting "samtools-index" { coord = _: true; } input);

stdenv.mkDerivation {
  name = "samtools-index";
  buildInputs = [ samtools ];
  buildCommand = ''
    ln -s ${input} input.bam
    samtools index -@ $NIX_BUILD_CORES ${optionalString (flags != null) flags} input.bam
    cp input.bam.bai $out
  '';
}
