{ bionix, nixpkgs }:

with bionix;

{
  align = callBionixE ./bwa-mem.nix;
  index = callBionixE ./bwa-index.nix;
}
