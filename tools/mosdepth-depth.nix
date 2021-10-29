{ bionix
, indexAttrs ? { }
, flags ? null
}:

with bionix;
with lib;

input:

stage {
  name = "mosdepth-depth";
  buildInputs = with pkgs; [ mosdepth ];
  buildCommand = ''
    mkdir $out
    ln -s ${input} input.bam
    ln -s ${bionix.samtools.index indexAttrs input} input.bam.bai
    mosdepth -t $NIX_BUILD_CORES ${optionalString (flags != null) flags} $out/out input.bam
  '';
  passthru.multicore = true;
}
