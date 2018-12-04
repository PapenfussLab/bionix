{ bionix, nixpkgs }:

with bionix;

{
  align = callBionix ./bowtie-align.nix;
  index = callBionix ./bowtie-index.nix;
}
