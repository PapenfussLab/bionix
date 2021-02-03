{ bionix }:

with bionix;

let
  gen = callBionixE ./sambamba-generic.nix;

in {
  /* Sort aligned reads
  Type: { nameSort :: bool, ... } -> bam -> bam
  */
  sort = callBionixE ./sambamba-sort.nix;

  /* Build an index
  Type: { ... } -> bam -> index
  */
  index = def gen {tool = "index"; };

  /* Merge bam files
  Type: { ... } -> [bam] -> bam
  */
  merge = def gen {tool = "merge"; };

  /* Slice a region out of a bam file
  Type: { region, ... } -> bam -> bam
  */
  slice = def gen {tool = "slice"; };

  /* Compute flag statistics
  Type: { ... } -> bam -> flagstat
  */
  flagstat = def gen {tool = "flagstat"; };

  /* Mark duplicates
  Type: { ... } -> bam -> bam
  */
  markdup = def gen {tool = "markdup"; };
}
