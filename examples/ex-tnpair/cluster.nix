{ bionix ? import <bionix> { }, tmpDir ? null }:

let
  bionix' = (bionix."${if tmpDir == null then "slurm" else "qsub"}" {
    ppn = 24;
    mem = 7;
    walltime = "3:00:00";
  } // bionix.lib.optionalAttrs (tmpDir != null) { inherit tmpDir; }).extend
    (self: super:
      with self; {
        minimap2.align = def super.minimap2.align {
          mem = 15;
          walltime = "16:00:00";
        };
        samtools = super.samtools // (with super.samtools; {
          markdup = def markdup { walltime = "12:00:00"; };
          fixmate = def fixmate { walltime = "10:00:00"; };
          sort = def sort {
            mem = 27;
            flags = "-m 1G";
          };
        });
      });

in import ./. { bionix = bionix'; }
