{ configuration, pkgs, lib ? pkgs.stdenv.lib }:

let
  pkgsModule = { config, lib, ... }: {
    config = { _module.args.pkgs = lib.mkDefault pkgs; };
  };

  modules = [ ./bwa.nix ./base.nix ];

  rawModules = lib.evalModules {
    modules = [ configuration ] ++ modules ++ [ pkgsModule ];
    specialArgs = { modulesPath = builtins.toString ./.; };
  };

in {
  inherit (rawModules) options config;
  bionix = rawModules.config.bionix;
}
