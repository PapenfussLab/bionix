{ config, lib, pkgs, ...}:

let cfg = config.samtools;
in with config.bionix;
with lib;

{
  options.samtools.sort = {
    mem = mkOption {
      type = types.str;
      default = "768M";
    };

    kmer = mkOption {
      type = types.int;
      default = 0;
    };

    pg = mkOption {
      type = types.bool;
      default = false;
    };

    order = mkOption {
      type = types.enum [ "coord" "name"];
      default = "coord";
    };
  };

  options.bionix.samtools.sort = mkOption {
    type = types.anything;
    description = "sorts an alignment";
  };

  config.bionix.samtools.sort = input: stage {
    name = "samtools-sort";
    buildInputs = with pkgs; [ samtools ];
    buildCommand = ''
      samtools sort -@ $NIX_BUILD_CORES \
        -m ${cfg.sort.mem} \
        ${optionalString (cfg.sort.kmer > 0) "-M -K ${toString cfg.sort.kmer}"} \
        ${optionalString (cfg.sort.order == "name") "-n"} \
        ${optionalString (cfg.sort.pg == false) "--no-PG"} \
        ${input} > $out
    '';
  };
}

