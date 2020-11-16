{ bionix
, nameSort ? false
, flags ? null
}:

input:

with bionix;
with lib;

let
  inherit (bionix.types) matchFiletype coordSort matchFileSorting;
in

assert (matchFiletype "sambamba-sort" { bam = _: true; } input);

stage {
  name = "sambamba-sort";
  buildInputs = [ pkgs.sambamba ];
  buildCommand = ''
    sambamba sort -t $NIX_BUILD_CORES \
      ${optionalString nameSort "-n"} \
      ${optionalString (flags != null) flags} \
      -o $out \
      ${input}
  '';
  passthru.filetype = if nameSort then bionix.types.nameSort input.filetype else coordSort input.filetype;
  passthru.multicore = true;
}
