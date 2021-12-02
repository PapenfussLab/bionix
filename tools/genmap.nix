{ bionix }:

with bionix;

{
  index = callBionixE ./genmap-index.nix;
  calcmap = callBionixE ./genmap-calcmap.nix;
}

