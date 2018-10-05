{ bionix
, nixpkgs
, flags ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-faidx" { fa = _: true; } input);

stdenv.mkDerivation {

  name = "samtools-faidx";
  buildInputs = [ samtools ];
  buildCommand = ''
    ln -s ${input} input.fasta
    samtools faidx ${optionalString (flags != null) flags} input.fasta
    cp input.fasta.fai $out
  '';
}
