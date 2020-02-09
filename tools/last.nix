{ bionix }:

with bionix;

rec {
  align = callBionixE ./last-align.nix;
  index = callBionixE ./last-index.nix;
}
