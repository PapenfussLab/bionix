{bionix}:

with bionix;

{
  app = pkgs.callPackage ./mutect-app.nix {inherit (pkgs) stdenv fetchurl makeWrapper unzip fetchFromGitHub;};
  call = callBionixE ./mutect-call.nix;
}
