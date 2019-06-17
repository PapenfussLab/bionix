{stdenv, lib, writeScript, coreutils}:

with lib;

{ ppn, mem, walltime, partition ? null, slurmFlags ? null, salloc ? "/usr/bin/salloc", srun ? "/usr/bin/srun" }:
drv:
  let ppnReified = if drv.multicore then ppn else 1;
  in lib.overrideDerivation drv ({ args, builder, name, ... }: {
    builder = stdenv.shell;
    args = let
      script = writeScript "slurm-script" ''
        #!${stdenv.shell}
        ${builder} ${lib.escapeShellArgs args}
      '';

      slurm = writeScript "slurm" ''
        #!${stdenv.shell}
        NIX_BUILD_CORES=${toString ppnReified}

        ${salloc} -c $NIX_BUILD_CORES --mem=${toString mem}G -t ${walltime} \
          -J "${name}" \
          ${optionalString (partition != null) "-p ${partition}"} \
          ${optionalString (slurmFlags != null) slurmFlags} \
          ${srun} ${script}
      '';

      in [ "-c" slurm ];
  })
