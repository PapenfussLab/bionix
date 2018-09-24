{ stdenv
, callPackage
, lib
, samtools
, flags ? null
}:

input:

with lib;

stdenv.mkDerivation {
  name = "samtools-faidx";
  buildInputs = [ samtools ];
  buildCommand = ''
    ln -s ${input} input.fasta
    samtools faidx ${optionalString (flags != null) flags} input.fasta
    cp input.fasta.fai $out
  '';
}