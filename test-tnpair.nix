with (import <nixpkgs> {});
with lib;

let
  bionix = (import <bionix> {}).extend (self: super: with self; {
    bwa = with super.bwa; {
      align = align;
      index = def index { flags = "-a is"; };
    };
  });

in

with bionix;

let
  fetchlocal = path: stdenv.mkDerivation {
    name = baseNameOf path;
    buildCommand = "ln -s ${path} $out";
  };
  fetchfq = attrs: types.tagFiletype (types.filetype.fq {}) (fetchlocal attrs);
  fetchfa = attrs: types.tagFiletype (types.filetype.fa {}) (fetchlocal attrs);

  ref = fetchfa ./example/ref.fa;

  alignWithRG = rg: bwa.align { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};
  sort = samtools.sort {};
  flagstat = samtools.flagstat {};
  check = fastqc.check {};
  callVariants = strelka.call {};

  tnpair = {
    tumour = {name = "mysample1"; files = {
        input1 = fetchfq ./example/sample1-1.fq;
        input2 = fetchfq ./example/sample1-2.fq;
      };
    };
    normal = {name = "mysample2"; files = {
        input1 = fetchfq ./example/sample2-1.fq;
        input2 = fetchfq ./example/sample2-2.fq;
      };
    };
  };

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
      ln -s ${bowtie.align {inherit ref;} tnpair.normal.files} $out/alignments/bowtie-normal.bam
      ln -s ${gridss.callVariants {} (with tnpairResult.alignments; [normal tumour])} $out/gridss
      ln -s ${gridss.call (with tnpairResult.alignments; [normal tumour])} $out/gridss2
      ln -s ${samtools.merge {} [tnpairResult.alignments.tumour tnpairResult.alignments.normal]} $out/alignments/merged.bam
      ln -s ${samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.tumour)} $out/alignments/${tnpair.tumour.name}.cram
      ln -s ${samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.normal)} $out/alignments/${tnpair.normal.name}.cram
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
