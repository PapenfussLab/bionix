{pkgs ? import <nixpkgs> {}}:

with pkgs;

stdenv.mkDerivation {
  name = "conda";
  src = ./conda.nix;
  phases = [ "buildPhase" ];
  buildInputs = [ conda ];
  buildPhase = ''
    HOME=$TMPDIR
    conda-shell-4.3.31 << EOF
    conda-install
    conda config --add channels bioconda
    conda install -y bwa
    bwa
    EOF
  '';
}
