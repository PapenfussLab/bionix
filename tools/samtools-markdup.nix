{ bionix
, nixpkgs
, flags ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-markdup" { bam = _: true; } input);
assert (matchFileSorting "samtools-markdup" { coord = _: true; } input);

stdenv.mkDerivation {
  name = "samtools-markdup";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools markdup ${optionalString (flags != null) flags} ${input} $out
  '';
  passthru.filetype = input.filetype;
}
