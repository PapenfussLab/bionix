{ bionix, nixpkgs }:

with nixpkgs;

{
  check = attrs: callPackage ./fastqc-check.nix attrs;
}
