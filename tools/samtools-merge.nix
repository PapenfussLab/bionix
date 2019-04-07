{ bionix
, flags ? null
, outfmt ? null
}:

inputs:

with bionix;
with lib;

let
  inherit (bionix.types) matchFiletype matchSorting;
  inputIsHomogenous = length (unique (map (matchFiletype "samtools-merge" {bam = x: x // {sorting = matchSorting "samtools-merge" {coord = _: "coord";} x;};}) inputs)) == 1;
in

assert inputIsHomogenous;

stage {
  name = "samtools-merge";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    samtools merge ${optionalString (flags != null) flags} out.bam ${concatStringsSep " " inputs}

    # Merge is non-deterministic with PG lines; if files have clashing PG IDs then a random
    # suffix is appended to make it unique. PG lines are stripped in the following to
    # resolve the issue.
    samtools reheader <(samtools view -H out.bam | grep -v '@PG') out.bam > $out
  '';
  passthru.filetype = (builtins.elemAt inputs 0).filetype;
}
