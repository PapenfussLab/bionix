{ config, lib, pkgs, ... }:

let cfg = config.bwa;

in with config.bionix;
with lib;

{
  options.bwa.align = {
    ref = mkOption {
      type = types.path;
      description = "the reference to align against";
    };

    bamOutput = mkOption {
      type = types.bool;
      default = true;
      description = "convert the output of BWA from SAM to BAM";
    };

    flags = mkOption {
      type = types.str;
      default = "";
      description = "additional flags to pass to BWA";
    };
  };

  options.bwa.index = {
    flags = mkOption {
      type = types.str;
      default = "";
      description = "additional flags to pass to BWA";
    };
  };

  options.bionix.bwa.index = mkOption {
    type = types.anything;
    description = "index a reference";
  };

  options.bionix.bwa.align = mkOption {
    type = types.anything;
    description = "align single or paird-end fastqs against a reference";
  };

  config.bionix.bwa.align = { input1, input2 ? null }:
    stage {
      name = "bwa-mem";
      buildInputs = with pkgs; [ bwa ] ++ optional cfg.align.bamOutput samtools;
      buildCommand = ''
        ln -s ${cfg.align.ref} ref.fa
        for f in ${bwa.index cfg.align.ref}/* ; do
          ln -s $f
        done
        bwa mem ${cfg.align.flags} -t $NIX_BUILD_CORES ref.fa ${input1} \
          ${optionalString (input2 != null) input2} \
          ${optionalString cfg.align.bamOutput "| samtools view -b"} \
          > $out
      '';
    };

  config.bionix.bwa.index = fa:
    stage {
      name = "bwa-index";
      buildInputs = with pkgs; [ bwa ];
      buildCommand = ''
        ln -s ${fa} ref.fa
        bwa index ${cfg.index.flags} ref.fa
        mkdir $out
        cp ref.fa.* $out
      '';
    };
}
