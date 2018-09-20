{ stdenv
, lib
, bwa
}:

ref:

with lib;

stdenv.mkDerivation {
  name = "bwa-index";
  buildInputs = [ bwa ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    bwa index ref.fa
    mkdir $out
    mv ref.fa.* $out
  '';
}
