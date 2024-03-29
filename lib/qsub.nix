{ stdenv, lib, writeScript }:

with lib;

let escape = x: if builtins.typeOf x == "string" then escapeShellArg x else x;

in
{ ppn
, mem
, walltime
, queue ? null
, qsubFlags ? null
, tmpDir
, sleepTime
, qsubPath ? "/usr/bin"
}:
drv:
let ppnReified = if drv.multicore then ppn else 1;
in
overrideDerivation drv ({ args, builder, name, ... }: {
  builder = "/bin/bash";
  args =
    let
      script = writeScript "qsub-script" ''
        #!${stdenv.shell}
        while [ ! -e ${tmpDir}/qsub-$PBS_JOBID ] ; do
          sleep ${toString sleepTime}
        done
        set -a
        . ${tmpDir}/qsub-$PBS_JOBID/nix-set
        set +a
        TMPDIR=${tmpDir}/qsub-$PBS_JOBID
        TEMP=$TMPDIR
        TMP=$TMPDIR
        NIX_BUILD_TOP=$TMPDIR
        cd $TMPDIR
        ${builder} ${concatMapStringsSep " " escape args} &> qsub-log
        echo $? > qsub-exit
      '';

      qsub = writeScript "qsub" ''
        #!${stdenv.shell}
        PATH=${qsubPath}
        SHELL=/bin/sh
        NIX_BUILD_CORES=${toString ppnReified}

        while : ; do
          qsub -l nodes=1:ppn=${toString ppnReified},mem=${
            toString mem
          }gb,walltime=${walltime} \
            -N "${name}" \
            ${optionalString (queue != null) "-q ${queue}"} \
            ${optionalString (qsubFlags != null) qsubFlags} \
            ${script} 2>&1 > id
          if [ $? -eq 0 ] ; then
            break
          fi
          if ! grep "Please retry" id > /dev/null ; then
            cat id >&2
            exit 1
          fi
          sleep ${toString sleepTime}
        done
        id=$(cat id)
        echo $id

        function cleanup {
          qdel $id 2>/dev/null || true
          sleep ${toString sleepTime}
          rm -rf ${tmpDir}/qsub-$id
        }
        trap cleanup INT TERM EXIT

        cp -r $TMPDIR ${tmpDir}/qsub-$id
        set > ${tmpDir}/qsub-$id/nix-set
        until qstat -f ''${id%%.} 2>&1 | grep "\(Unknown Job\|job_state = C\)" > /dev/null ; do
          sleep ${toString sleepTime}
        done
        cat ${tmpDir}/qsub-$id/qsub-log
        if [ -e ${tmpDir}/qsub-$id/qsub-exit ]; then
          exitCode=$(cat ${tmpDir}/qsub-$id/qsub-exit)
        else
          exitCode=1
        fi
        exit $exitCode
      '';

    in
    [ "-c" qsub ];
})
