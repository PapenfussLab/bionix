{ bionix }:

with bionix;

rec {
  /* Align read against a reference: defaults to bwa-mem */
  align = bwa-mem;

  /* Align reads against a reference using bwa-mem
  Type: bwa-mem :: {ref :: fasta, bamOutput :: bool, ...} -> {input1, input2} -> bam/sam
  */
  bwa-mem = callBionixE ./bwa-mem.nix;
  /* Creates an reference index for BWA
  Type: index :: {...} -> fasta -> BWA index
  */
  index = callBionixE ./bwa-index.nix;
}
