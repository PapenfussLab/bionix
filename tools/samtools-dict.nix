{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-dict" { fa = _: true; } input);

stage {
  name = "samtools-dict";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    samtools dict ${optionalString (flags != null) flags} ${input} > $out
  '';
}
