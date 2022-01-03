{ bionix
}:

input:

with bionix;
with lib;

stage {
  name = "samtools-flagstat";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = "samtools flagstat -@ $NIX_BUILD_CORES ${input} > $out";
  passthru.multicore = true;
}
