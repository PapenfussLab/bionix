{lib}:

with lib;


{ ppn ? 1, mem ? 1, walltime ? "24:00:00" }: drv: lib.overrideDerivation drv ({ args, builder, ... }: {
  builder = "/bin/bash";
  args = let
    script = writeScript "qsub-script" ''
      #!${stdenv.shell}
      while [ ! -e /stornext/HPCScratch/$PBS_JOBID ] ; do
        sleep 5
      done
      set -a
      . /stornext/HPCScratch/$PBS_JOBID
      set +a
      TMPDIR=/tmp/$PBS_JOBID
      TEMP=$TMPDIR
      TMP=$TMPDIR
      NIX_BUILD_TOP=$TMPDIR
      mkdir $TMPDIR
      cd $TMPDIR
      rm /stornext/HPCScratch/$PBS_JOBID
      ${builder} ${lib.escapeShellArgs args} && touch /stornext/HPCScratch/$PBS_JOBID
      cd /
      rm -rf $TMPDIR
    '';

    qsub = writeScript "qsub" ''
      #!/bin/bash
      PATH=/usr/bin:/bin:/usr/sbin:/sbin
      SHELL=/bin/sh
      NIX_BUILD_CORES=${toString ppn}
      id=$(qsub -l nodes=1:ppn=${toString ppn},mem=${toString mem}gb,walltime=${walltime} ${script})
      set > /stornext/HPCScratch/$id
      while qstat ''${id%%.} 2> /dev/null > /dev/null ; do
        sleep 5
      done
      if [[ -e /stornext/HPCScratch/$id ]] ; then
        rm /stornext/HPCScratch/$id
        exit 0
      fi
      exit 1
    '';

    in [ "-c" qsub ];
})
