{bionix}:

with bionix;

input:

stage {
  name = "gzip";
  buildInputs = with pkgs; [ gzip ];
  buildCommand = "gzip < ${input} > $out";
  passthru.filetype = filetype.gzip input.filetype;
}
