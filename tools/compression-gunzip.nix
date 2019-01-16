{bionix}:

with bionix;

input:

stage {
  name = "gunzip";
  buildInputs = with pkgs; [ gzip ];
  buildCommand = "gunzip < ${input} > $out";
  passthru.filetype = types.gunzip input;
}
