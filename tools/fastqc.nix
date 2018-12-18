{ bionix, nixpkgs }:

with nixpkgs;
with bionix;

{
  fastqc = callPackage ./fastqc-app.nix {};
  check = callBionixE ./fastqc-check.nix;
}
