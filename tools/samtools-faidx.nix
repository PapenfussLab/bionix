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
    samtools faidx ${optionalString (flags != null) flags} ${input} > $out
  '';
}
