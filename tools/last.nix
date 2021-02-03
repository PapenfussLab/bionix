{ bionix }:

with bionix;

rec {
  /* Align reads against a reference genome
  Type: { ref :: fastq, ... } -> { input1, input2} -> last
  */
  align = callBionixE ./last-align.nix;
  index = callBionixE ./last-index.nix;
}
