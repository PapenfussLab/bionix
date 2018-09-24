{pkgs ? import <nixpkgs> {}}:

with pkgs;
with lib;

let
  ref = ./example/ref.fa;
  alignWithRG = rg: callPackage ./tools/bwa.nix { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};
  sort = callPackage ./tools/samtools-sort.nix { };
  callVariants = callPackage ./tools/strelka.nix { inherit ref; };

  tnpair = { tumour = {name = "mysample1"; files = {input1 = ./example/sample1-1.fq; input2 = ./example/sample1-2.fq;};};
               normal = {name = "mysample2"; files = {input1 = ./example/sample2-1.fq; input2 = ./example/sample2-1.fq;};};};

  processPair = { tumour, normal }: rec {
    alignments = { normal = sort(alignWithRG normal.name normal.files); tumour = sort (alignWithRG tumour.name tumour.files); };
    variants = callVariants alignments;
  };

  #results = map processPair tnpairs;
  tnpairResult = processPair tnpair;

  testNaming = stdenv.mkDerivation {
    name = "test-naming";
    buildCommand = ''
      mkdir $out
      ln -s ${tnpairResult.variants} $out/strelka
      mkdir $out/alignments
      ln -s ${tnpairResult.alignments.tumour} $out/alignments/${tnpair.tumour.name}.bam
      ln -s ${tnpairResult.alignments.normal} $out/alignments/${tnpair.normal.name}.bam
    '';
  };

in testNaming
