{ bionix }:

with bionix;

{
  view = callBionixE ./samtools-view.nix;
  faidx = callBionixE ./samtools-faidx.nix;
  flagstat = callBionixE ./samtools-flagstat.nix;
  index = callBionixE ./samtools-index.nix;
  dict = callBionixE ./samtools-dict.nix;
  sort = callBionixE ./samtools-sort.nix;
  merge = callBionixE ./samtools-merge.nix;
  markdup = callBionixE ./samtools-markdup.nix;
  fixmate = callBionixE ./samtools-fixmate.nix;
}
