with (import <nixpkgs> {});
with (import <bionix> {});
with lib;

let
  ref = ./example/ref.fa;
  alignWithRG = rg: bwa.align { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};
  sort = samtools.sort { };
  flagstat = samtools.flagstat {};
  check = fastqc.check {};
  callVariants = strelka.call { inherit ref; };

  tnpair = { tumour = {name = "mysample1"; files = {input1 = ./example/sample1-1.fq; input2 = ./example/sample1-2.fq;};};
             normal = {name = "mysample2"; files = {input1 = ./example/sample2-1.fq; input2 = ./example/sample2-1.fq;};};};

  processPair = { tumour, normal }: rec {
    alignments = mapAttrs (_: x: sort (alignWithRG x.name x.files)) { inherit normal tumour; };
    variants = callVariants alignments;
  };

  tnpairResult = processPair tnpair;

  testNaming = stdenv.mkDerivation {
    name = "test-naming";
    buildCommand = ''
      mkdir $out
      ln -s ${tnpairResult.variants} $out/strelka
      mkdir $out/alignments
      ln -s ${tnpairResult.alignments.tumour} $out/alignments/${tnpair.tumour.name}.bam
      ln -s ${tnpairResult.alignments.normal} $out/alignments/${tnpair.normal.name}.bam
      ln -s ${flagstat tnpairResult.alignments.tumour} $out/alignments/${tnpair.tumour.name}.flagstat
      ln -s ${flagstat tnpairResult.alignments.normal} $out/alignments/${tnpair.normal.name}.flagstat
      mkdir $out/fastqc
      ln -s ${check tnpair.tumour.files.input1} $out/fastqc/${tnpair.tumour.name}.1
      ln -s ${check tnpair.tumour.files.input2} $out/fastqc/${tnpair.tumour.name}.2
      ln -s ${check tnpair.normal.files.input1} $out/fastqc/${tnpair.normal.name}.1
      ln -s ${check tnpair.normal.files.input2} $out/fastqc/${tnpair.normal.name}.2
    '';
  };

in testNaming
