{ tmpDir ? null }:

let
  bionix = import <bionix> {
    overlays = [
      (_: super:
        super."${if tmpDir == null then "slurm" else "qsub"}"
          {
            ppn = 24;
            mem = 7;
            walltime = "3:00:00";
          } // super.lib.optionalAttrs (tmpDir != null) { inherit tmpDir; })

      (self: super:
        with super; {
          minimap2.align = def minimap2.align {
            mem = 15;
            walltime = "16:00:00";
          };
          samtools = with samtools; {
            markdup = def markdup { walltime = "12:00:00"; };
            fixmate = def fixmate { walltime = "10:00:00"; };
            sort = def sort {
              mem = 27;
              flags = "-m 1G";
            };
          };
        })
    ];
  };
in
import ./. { inherit bionix; }
