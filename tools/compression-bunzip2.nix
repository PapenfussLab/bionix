{bionix}:

with bionix;

input:

stage {
  name = "bunzip2";
  buildInputs = with pkgs; [ bzip2 ];
  buildCommand = "bunzip2 < ${input} > $out";
  passthru.filetype = types.bunzip2 input;
}
