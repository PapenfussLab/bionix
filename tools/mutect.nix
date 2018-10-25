{bionix, nixpkgs}:

with nixpkgs;
with bionix;

{
  app = callPackage ./mutect-app.nix {inherit (nixpkgs) stdenv fetchurl makeWrapper unzip fetchFromGitHub;};
  call = callBionix ./mutect-call.nix;
}
