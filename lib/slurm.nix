{ stdenv, lib, writeScript, coreutils }:

with lib;

let escape = x: if builtins.typeOf x == "string" then escapeShellArg x else x;

in
{ ppn
, mem
, walltime
, partition ? null
, slurmFlags ? null
, salloc ? "/usr/bin/salloc"
, srun ? "/usr/bin/srun"
, ...
}:
drv:
let ppnReified = if drv.multicore then ppn else 1;
in
overrideDerivation drv ({ args, builder, name, ... }: {
  builder = stdenv.shell;
  args =
    let
      slurm = writeScript "slurm" ''
        #!${stdenv.shell}
        NIX_BUILD_CORES=${toString ppnReified}

        ${salloc} -c $NIX_BUILD_CORES --mem=${toString mem}G -t ${walltime} \
          -J "${name}" \
          ${optionalString (partition != null) "-p ${partition}"} \
          ${optionalString (slurmFlags != null) slurmFlags} \
          ${srun} ${builder} ${concatMapStringsSep " " escape args}
      '';

    in
    [ "-c" slurm ];
})
