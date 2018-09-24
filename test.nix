{pkgs ? import <nixpkgs> {}}:

with pkgs;
with lib;

let
  ref = ./example/ref.fa;
  alignWithRG = rg: callPackage ./tools/bwa.nix { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};
  sort = callPackage ./tools/samtools-sort.nix { };
  callVariants = callPackage ./tools/platypus.nix { inherit ref; };

  samples = [ {name = "mysample1"; files = {input1 = ./example/sample1-1.fq; input2 = ./example/sample1-2.fq;};}
             {name = "mysample2"; files = {input1 = ./example/sample2-1.fq; input2 = ./example/sample2-1.fq;};} ];

  alignments = map (i: sort (alignWithRG i.name i.files)) samples;
  variants = callVariants alignments;

  testNaming = stdenv.mkDerivation {
    name = "test-naming";
    buildCommand = ''
      mkdir $out
      ln -s ${variants} $out/myfancyname
      mkdir $out/alignments
      ${concatStringsSep "\n" (zipListsWith (s: a: "ln -s ${a} $out/alignments/${s.name}.bam") samples alignments)}
    '';
  };

in testNaming
