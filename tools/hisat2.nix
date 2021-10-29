{ bionix }:

with bionix;

rec {
  /* Align read against a reference
    Type: hisat2-mem :: {ref :: fasta, bamOutput :: bool, ...} -> {input1, input2} -> bam/sam
  */
  align = callBionixE ./hisat2-align.nix;

  /* Creates an reference index for HISAT2
    Type: index :: {...} -> fasta -> HISAT2 index
  */
  index = callBionixE ./hisat2-index.nix;

  extractSpliceSites = callBionixE ./hisat2-extractSpliceSites.nix;
  extractExons = callBionixE ./hisat2-extractExons.nix;
}
