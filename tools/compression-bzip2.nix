{ bionix }:

with bionix;

input:

stage {
  name = "bzip2";
  buildInputs = with pkgs; [ bzip2 ];
  buildCommand = "bzip2 < ${input} > $out";
  passthru.filetype = filetype.bz2 input.filetype;
}
