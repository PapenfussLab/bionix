{ bionix, nixpkgs }:

with bionix;

{
  call = callBionixE ./platypus-callVariants.nix;
}
