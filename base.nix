{ config, pkgs, lib, ... }:

with lib;

let
  strip = drv:
    let
      stripCommand = ''

        function rewrite {
          sed -i 's|[A-Za-z0-9+/]\{32\}-bionix|00000000000000000000000000000000-bionix|g' $1
        }
        function rewriteOutput {
          if [ -f ''${!1} ] ; then
            rewrite ''${!1}
          else
            for f in $(find ''${!1} -type f) ; do
              rewrite $f
            done
          fi
        }
        for o in $outputs ; do
          rewriteOutput $o
        done
      '';
    in drv.overrideAttrs (attrs:
      if attrs ? buildCommand then {
        buildCommand = attrs.buildCommand + stripCommand;
      } else {
        fixupPhase = (if attrs ? fixupPhase then attrs.fixupPhase else "")
          + stripCommand;
      });

in {
  options.bionix.stage = mkOption {
    internal = true;
    type = types.anything;
    description = "stage";
  };

  config.bionix.stage =
    x@{ name, stripStorePaths ? true, multicore ? false, ... }:
    (if stripStorePaths then strip else x: x) (pkgs.stdenvNoCC.mkDerivation (x
      // {
        name = "bionix-" + name;
        inherit multicore;
      }));
}
