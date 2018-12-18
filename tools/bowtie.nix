{ bionix, nixpkgs }:

with bionix;

{
  align = callBionixE ./bowtie-align.nix;
  index = callBionixE ./bowtie-index.nix;
}
