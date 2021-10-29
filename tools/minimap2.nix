{ bionix }:

with bionix;

rec {
  /* Align read against a reference 
    * Type: align :: {ref :: fasta, bamOutput :: bool, preset :: string, ...} -> {input1, input2} -> bam/sam
  */
  align = callBionixE ./minimap2-align.nix;
}
