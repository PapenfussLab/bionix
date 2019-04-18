{stdenv, lib, writeScript, coreutils}:

with lib;

{ ppn, mem, walltime, partition ? null, slurmFlags ? null}:
drv:
  let ppnReified = if drv.multicore then ppn else 1;
  in lib.overrideDerivation drv ({ args, builder, name, ... }: {
    builder = stdenv.shell;
    args = let
      script = writeScript "sbatch-script" ''
        #!${stdenv.shell}
        ${builder} ${lib.escapeShellArgs args}
      '';

      sbatch = writeScript "sbatch" ''
        #!${stdenv.shell}
        NIX_BUILD_CORES=${toString ppnReified}

        salloc -c $NIX_BUILD_CORES --mem=${toString mem}G -t ${walltime} \
          -J "${name}" \
          ${optionalString (partition != null) "-p ${partition}"} \
          ${optionalString (slurmFlags != null) slurmFlags} \
          ${script}
      '';

      in [ "-c" sbatch ];
  })
