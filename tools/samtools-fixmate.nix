{ bionix
, mateScore ? true
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-fixmate" { bam = _: true; } input);
assert (matchFileSorting "samtools-fixmate" { name = _: true; } input);

stage {
  name = "samtools-fixmate";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    samtools fixmate \
      ${optionalString mateScore "-m"} \
      ${optionalString (flags != null) flags} -O bam ${input} $out
  '';
  passthru.filetype = input.filetype;
}
