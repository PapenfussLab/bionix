{ bionix
, k ? 30
, e ? 2
}:

with bionix;

ref: stage {
  name = "genmap-mappability";
  buildInputs = with pkgs; [ genmap ];
  buildCommand = ''
    mkdir $out
    genmap map \
      -K ${toString k} \
      -E ${toString e} \
      --wig \
      -I ${bionix.genmap.index {} ref} \
      -O $out \
      -T $NIX_BUILD_CORES
  '';
  passthru.multicore = true;
}
