{ bionix
, nixpkgs
, mateScore ? false
, flags ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-fixmate" { bam = _: true; } input);
assert (matchFileSorting "samtools-fixmate" { name = _: true; } input);

stdenv.mkDerivation {
  name = "samtools-fixmate";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools fixmate \
      ${optionalString mateScore "-m"} \
      ${optionalString (flags != null) flags} -O bam ${input} $out
  '';
  passthru.filetype = input.filetype;
}
