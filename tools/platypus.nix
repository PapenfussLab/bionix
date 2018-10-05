{ bionix, nixpkgs }:

with bionix;

{
  call = callBionix ./platypus-callVariants.nix;
}
