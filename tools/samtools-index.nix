{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-index" { bam = _: true; } input);
assert (matchFileSorting "samtools-index" { coord = _: true; } input);

stage {
  name = "samtools-index";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    ln -s ${input} input.bam
    samtools index -@ $NIX_BUILD_CORES ${optionalString (flags != null) flags} input.bam
    cp input.bam.bai $out
  '';
}
