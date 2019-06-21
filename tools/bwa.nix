{ bionix }:

with bionix;

rec {
  app2 = pkgs.callPackage ./bwa-mem2-app.nix {};

  /* Align read against a reference: defaults to bwa-mem */
  align = bionix.bwa.mem;

  /* Align reads against a reference using bwa-mem
  Type: bwa-mem :: {ref :: fasta, bamOutput :: bool, ...} -> {input1, input2} -> bam/sam
  */
  mem = callBionixE ./bwa-mem.nix;
  mem2 = callBionixE ./bwa-mem2.nix;
  /* Creates an reference index for BWA
  Type: index :: {...} -> fasta -> BWA index
  */
  index = callBionixE ./bwa-index.nix;
  index2 = callBionixE ./bwa-index2.nix;
}
