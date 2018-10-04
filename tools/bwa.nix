{ bionix, nixpkgs }:

with bionix;

{
  align = callBionix ./bwa-mem.nix;
  index = callBionix ./bwa-index.nix;
}
