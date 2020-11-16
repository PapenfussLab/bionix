{ bionix
, flags ? null
, tool
, region ? null
}:

input:

with bionix;
with lib;

let
  inherit (bionix.types) matchFiletype coordSort matchFileSorting;
in

assert (matchFiletype "sambamba-${tool}" { bam = _: true; } input);

stage {
  name = "sambamba-${tool}";
  buildInputs = [ pkgs.sambamba ];
  buildCommand = ''
    sambamba ${tool} -t $NIX_BUILD_CORES \
      ${optionalString (flags != null) flags} \
      ${if tool == "merge" then "$out ${concatStringsSep " " input}" else if tool == "slice" then "${input} ${region} > $out" else if tool == "flagstat" then "${input} > $out" else "${input} $out"}
  '';
  passthru.filetype = if tool == "flagstat" || tool == "index" then null else if tool == "merge" then (head input).filetype else input.filetype;
  passthru.multicore = true;
}
