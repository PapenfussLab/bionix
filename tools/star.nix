{ bionix }:

with bionix;

{
  /* Align RNA against a reference genome
  Type: { ref, ... } -> { input1, input2 } -> bam
  */
  align = callBionixE ./star-align.nix;
  index = callBionixE ./star-index.nix;
}
