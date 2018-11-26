{stdenv, lib, writeScript}:

{ ppn ? 1, mem ? 1, walltime ? "24:00:00", tmpDir ? "/tmp" , sleepTime ? 60}: drv: lib.overrideDerivation drv ({ args, builder, name, ... }: {
  builder = "/bin/bash";
  args = let
    script = writeScript "qsub-script" ''
      #!${stdenv.shell}
      while [ ! -e ${tmpDir}/$PBS_JOBID ] ; do
        sleep ${toString sleepTime}
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
      #!${stdenv.shell}
      PATH=/usr/bin:/bin:/usr/sbin:/sbin
      SHELL=/bin/sh
      NIX_BUILD_CORES=${toString ppn}

      while : ; do
        qsub -l nodes=1:ppn=${toString ppn},mem=${toString mem}gb,walltime=${walltime} -N "${name}" ${script} 2>&1 > id
        if [ $? -eq 0 ] ; then
          break
        fi
        if ! grep "Please retry" id > /dev/null ; then
          exit 1
        fi
        sleep ${toString sleepTime}
      done
      id=$(cat id)

      function cleanup {
        qdel $id 2>/dev/null || true
        sleep ${toString sleepTime}
        rm -rf ${tmpDir}/$id
      }
      trap cleanup INT TERM EXIT

      cp -r $TMPDIR ${tmpDir}/$id
      set > ${tmpDir}/$id/nix-set
      until qstat -f ''${id%%.} 2>&1 | grep "\(Unknown Job\|job_state = C\)" > /dev/null ; do
        sleep ${toString sleepTime}
      done
      cat ${tmpDir}/$id/qsub-stderr >&2
      cat ${tmpDir}/$id/qsub-stdout
      if [ -e ${tmpDir}/$id/qsub-exit ]; then
        exitCode=$(cat ${tmpDir}/$id/qsub-exit)
      else
        exitCode=1
      fi
      exit $exitCode
    '';

    in [ "-c" qsub ];
})
