{ bionix }:

with bionix;

rec {
  app = pkgs.callPackage ./snap-app.nix {};

  /* Align reads against a reference
  Type: snap :: {ref :: fasta, bamOutput :: bool, ...} -> {input1, input2} -> bam/sam
  */
  align = callBionixE ./snap-align.nix;

  /* Creates an reference index for SNAP
  Type: index :: {...} -> fasta -> SNAP index
  */
  index = callBionixE ./snap-index.nix;
}
