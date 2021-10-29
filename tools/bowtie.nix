{ bionix }:

with bionix;

{
  /* Align a sequence against a reference
    Type: align :: {ref :: fasta, bamOutput :: bool, ...} -> {input1 :: fastq, input2 :: fastq} -> bam/sam
  */
  align = callBionixE ./bowtie-align.nix;

  /* Create a Bowtie index
    Type: index :: {seed :: int, ...} -> fasta -> index
  */
  index = callBionixE ./bowtie-index.nix;
}
