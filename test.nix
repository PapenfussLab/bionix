{pkgs ? import <nixpkgs> {}}:

with pkgs;

let

  qsub = drv: lib.overrideDerivation drv ({ ppn ? 1, mem ? 1, walltime ? "24:00:00", args, builder, ... }: {
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
  });

  qstat = stdenv.mkDerivation rec {
    name = "qstat";
    buildCommand = "/usr/bin/qstat  > $out";
  };

  grep = x: stdenv.mkDerivation rec {
    name = "grep";
    buildInputs = [ bwa ];
    buildCommand = "/usr/bin/grep distefano.l ${x} > $out";
  };

  dummy = x: stdenv.mkDerivation rec {
    name = "dummy";
    buildCommand = "sleep 5 && echo ${toString x} > $out";
  };


in stdenv.mkDerivation rec {
  name = "dummies";
  dummies = map dummy [1 2 3 4 5];
  buildCommand = "for f in ${toString dummies} ; cat $f >> $out";
}

