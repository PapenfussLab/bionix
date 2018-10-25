{ bionix
, nixpkgs
, flags ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-dict" { fa = _: true; } input);

stdenv.mkDerivation {
  name = "samtools-dict";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools dict ${optionalString (flags != null) flags} ${input} > $out
  '';
}
