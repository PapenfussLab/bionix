{stdenv, lib, writeScript}:

{ ppn ? 1, mem ? 1, walltime ? "24:00:00", tmpDir ? "/tmp" }: drv: lib.overrideDerivation drv ({ args, builder, ... }: {
  builder = "/bin/bash";
  args = let
    script = writeScript "qsub-script" ''
      #!${stdenv.shell}
      while [ ! -e ${tmpDir}/$PBS_JOBID ] ; do
        sleep 5
      done
      TMPDIR=${tmpDir}/$PBS_JOBID
      TEMP=$TMPDIR
      TMP=$TMPDIR
      NIX_BUILD_TOP=$TMPDIR
      cd $TMPDIR
      set -a
      . nix-set
      set +a
      ${builder} ${lib.escapeShellArgs args} > qsub-stdout 2> qsub-stderr
      echo $? > qsub-exit
    '';

    qsub = writeScript "qsub" ''
      #!/bin/bash
      PATH=/usr/bin:/bin:/usr/sbin:/sbin
      SHELL=/bin/sh
      NIX_BUILD_CORES=${toString ppn}
      id=$(qsub -l nodes=1:ppn=${toString ppn},mem=${toString mem}gb,walltime=${walltime} ${script})

      function cleanup {
        qstat ''${id%%.} 2> /dev/null > /dev/null && qdel $id || true
        sleep 5
        rm -rf ${tmpDir}/$id
      }
      trap cleanup INT TERM EXIT

      cp -r $TMPDIR ${tmpDir}/$id
      set > ${tmpDir}/$id/nix-set
      while qstat ''${id%%.} 2> /dev/null > /dev/null ; do
        sleep 5
      done
      cat ${tmpDir}/$id/qsub-stderr >&2
      cat ${tmpDir}/$id/qsub-stdout
      exitCode=$(cat ${tmpDir}/$id/qsub-exit)
      exit $exitCode
    '';

    in [ "-c" qsub ];
})
