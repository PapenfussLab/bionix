# This is a translation of the Nextflow example found at
# https://www.nextflow.io/example1.html
{ bionix ? import ./../.. {}
, input ? ./sample.fa}:

with bionix;
with lib;

let
  splitSequences = fa: stage {
    name = "splitSequences";
    buildInputs = [ pkgs.gawk ];
    buildCommand = ''
      awk '/^>/{f="seq_"++d} {print > f}' ${fa}
      mkdir $out
      cp seq* $out
    '';
  };

  reverse = fa: stage {
    name = "reverse";
    buildCommand = ''
      ${pkgs.utillinux}/bin/rev ${fa} > $out
    '';
  };

in pipe [
  splitSequences
  (each reverse)
] input
