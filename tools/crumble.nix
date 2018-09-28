{ bionix, nixpkgs }:

with nixpkgs;

{
  toCram = attrs: callPackage ./crumble-toCram.nix attrs;
}
