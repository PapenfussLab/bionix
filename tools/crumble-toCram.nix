{ stdenv
, lib
, callPackage
, crumble ? callPackage ./crumble-app.nix {}
, flags ? null}:

with lib;

input:

stdenv.mkDerivation {
  name = "crumble";
  buildInputs = [ crumble ];
  buildCommand = "crumble ${optionalString (flags != null) flags} ${input} $out";
}
