{ bionix }:

with bionix;

{
  align = callBionixE ./star-align.nix;
  index = callBionixE ./star-index.nix;
}
