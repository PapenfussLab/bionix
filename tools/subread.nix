{ bionix }:

with bionix;

{
  align = callBionixE ./subread-align.nix;
  index = callBionixE ./subread-index.nix;
}
