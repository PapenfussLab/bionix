{ bionix
, nixpkgs
, flags ? null
, outfmt ? null
}:

inputs:

with nixpkgs;
with lib;

let
  inherit (bionix.types) matchFiletype option-sort;
  inputIsSorted = input: matchFiletype "samtools-merge" {
       bam = _: true; #{sorting, ...}: sorting == option-sort.some (bionix.types.sorting.coord {});
     } input;
in

assert (all inputIsSorted inputs);

stdenv.mkDerivation {
  name = "samtools-merge";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools merge ${optionalString (flags != null) flags} $out ${concatStringsSep " " inputs}
  '';
  passthru.filetype = (builtins.elemAt inputs 0).filetype;
}
