{ bionix
, nixpkgs
, flags ? null
}:

with nixpkgs;
with lib;

input:

stdenv.mkDerivation {
  name = "crumble";
  buildInputs = [ bionix.crumble.crumble ];
  buildCommand = "crumble ${optionalString (flags != null) flags} ${input} $out";
}
