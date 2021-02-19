{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-index" { bam = _: true; cram = _: true; } input);
assert (matchFileSorting "samtools-index" { coord = _: true; } input);

let
  ext = matchFiletype "samtools-index-ext" { bam = _: "bam"; cram = _: "cram"; };
in

stage {
  name = "samtools-index";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    ln -s ${input} input.${ext input}
    samtools index -@ $NIX_BUILD_CORES ${optionalString (flags != null) flags} input.${ext input} $out
  '';
  passthru.multicore = true;
}
