{ bionix, nixpkgs }:

with bionix;

{
  faidx = callBionix ./samtools-faidx.nix;
  flagstat = callBionix ./samtools-flagstat.nix;
  index = callBionix ./samtools-index.nix;
  sort = callBionix ./samtools-sort.nix;
}
