{ bionix }:

with bionix;

rec {
  index = callBionixE ./whisper-index.nix;

  /* Align reads against a reference
  Type: { ref, ... } -> { input1, input2 } -> bam
  */
  align = callBionixE ./whisper-align.nix;
}
