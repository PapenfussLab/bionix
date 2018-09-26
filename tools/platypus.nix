{ bionix, nixpkgs }:

with nixpkgs;

{
  call = attrs: callPackage ./platypus-callVariants.nix attrs;
}
