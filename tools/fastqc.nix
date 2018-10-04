{ bionix, nixpkgs }:

with nixpkgs;
with bionix;

{
  fastqc = callPackage ./fastqc-app.nix {};
  check = callBionix ./fastqc-check.nix;
}
