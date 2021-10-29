{ bionix
, flags ? null
, outfmt ? null
}:

inputs:

with bionix;
with lib;

let
  inherit (bionix.types) matchFiletype matchSorting;
  inputIsHomogenous = length (unique (map (matchFiletype "samtools-merge" { bam = x: x // { sorting = matchSorting "samtools-merge" { coord = _: "coord"; name = _: "name"; } x; }; }) inputs)) == 1;
  nameSorted = matchFiletype "samtools-merge" { bam = matchSorting "samtools-merge" { coord = _: false; name = _: true; }; } (lib.head inputs);
in

assert inputIsHomogenous;

if length inputs == 1 then head inputs else

stage {
  name = "samtools-merge";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    samtools merge ${optionalString (flags != null) flags} \
      ${if nameSorted then "-n" else ""} \
      out.bam ${concatStringsSep " " inputs}

    # Merge is non-deterministic with PG lines; if files have clashing PG IDs then a random
    # suffix is appended to make it unique. PG lines are stripped in the following to
    # resolve the issue.
    samtools reheader <(samtools view -H out.bam | grep -v '@PG') out.bam > $out
  '';
  passthru.filetype = (builtins.elemAt inputs 0).filetype;
}
