{ bionix
}:

input:

with bionix;
with lib;

stage {
  name = "samtools-index";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = "samtools flagstat -@ $NIX_BUILD_CORES ${input} > $out";
}
