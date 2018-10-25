{ bionix, nixpkgs }:

with bionix;

{
  view = callBionix ./samtools-view.nix;
  faidx = callBionix ./samtools-faidx.nix;
  flagstat = callBionix ./samtools-flagstat.nix;
  index = callBionix ./samtools-index.nix;
  dict = callBionix ./samtools-dict.nix;
  sort = callBionix ./samtools-sort.nix;
  merge = callBionix ./samtools-merge.nix;
}
