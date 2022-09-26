{
  bionix,
  indexAttrs ? {},
  targets ? null,
  flags ? null,
}:
with bionix;
with lib;
  input: let
    handleTarget = x: let
      type = builtins.typeOf x;
      handler = handlers."${type}" or (builtins.throw "mosdepth-depth:unhandled target type:${type}");
      handlers = {
        string = handleTarget [x];
        list = let
          file = pkgs.writeText "target.bed" (concatStringsSep "\n" x);
        in "-b ${file}";
        path = "-b ${x}";
        set = "-b ${x}";
      };
    in
      handler;
  in
    stage {
      name = "mosdepth-depth";
      buildInputs = with pkgs; [mosdepth];
      buildCommand = ''
        mkdir $out
        ln -s ${input} input.bam
        ln -s ${bionix.samtools.index indexAttrs input} input.bam.bai
        mosdepth -t $NIX_BUILD_CORES \
          ${optionalString (targets != null) (handleTarget targets)} \
          ${optionalString (flags != null) flags} $out/out input.bam
      '';
      passthru.multicore = true;
    }
