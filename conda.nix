{ stdenv, lib, writeScript, conda }:

{
  buildCondaEnv = { depends ? [], run }: stdenv.mkDerivation {
    name = "conda-env";
    buildInputs = [ conda ] ++ depends;
    buildCommand = ''
      mkdir $out
      HOME=$out
      conda-shell-4.3.31 << EOF
      conda-install
      ${run}
      EOF
    '';
  };

  withCondaEnv = env: run: stdenv.mkDerivation {
    name = "with-conda-env";
    buildCommand = ''
      #!${stdenv.shell}
      export HOME=${env}
      ${conda}/bin/conda-shell-4.3.31 << EOF
      ${run}
      EOF
    '';
  };
}

