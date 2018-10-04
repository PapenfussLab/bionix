{ bionix
, nixpkgs
, flags ? null
}:

with nixpkgs;
with lib;

input:

stdenv.mkDerivation {
  name = "fastqc-check";
  buildInputs = [ bionix.fastqc.fastqc ];
  buildCommand = ''
    mkdir $out
    fastqc \
      -o $out \
      ${optionalString (flags != null) flags} \
      ${input}
  '';
}
