{ bionix }:

with bionix;

{
  /* Mark duplicates
    Type: { ... } -> bam -> bam
  */
  markDuplicates = callBionixE ./picard-markDuplicates.nix;
}
