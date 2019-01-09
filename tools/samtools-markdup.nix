{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-markdup" { bam = _: true; } input);
assert (matchFileSorting "samtools-markdup" { coord = _: true; } input);

stage {
  name = "samtools-markdup";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    samtools markdup ${optionalString (flags != null) flags} ${input} $out
  '';
  passthru.filetype = input.filetype;
}
