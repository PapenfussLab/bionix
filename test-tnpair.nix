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

  testNaming = linkDrv [
    (ln tnpairResult.variants "strelka")
    (ln (bowtie.align {inherit ref;} tnpair.normal.files) "alignments/bowtie-normal.bam")
    (ln (gridss.callVariants {} (with tnpairResult.alignments; [normal tumour])) "gridss")
    (ln (gridss.call (with tnpairResult.alignments; [normal tumour])) "gridss2")
    (ln (samtools.merge {} [tnpairResult.alignments.tumour tnpairResult.alignments.normal]) "alignments/merged.bam")
    (ln (samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.tumour)) "alignments/${tnpair.tumour.name}.cram")
    (ln (samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.normal)) "alignments/${tnpair.normal.name}.cram")
    (ln (flagstat tnpairResult.alignments.tumour) "alignments/${tnpair.tumour.name}.flagstat")
    (ln (flagstat tnpairResult.alignments.normal) "alignments/${tnpair.normal.name}.flagstat")
    (ln (check tnpair.tumour.files.input1) "fastqc/${tnpair.tumour.name}.1")
    (ln (check tnpair.tumour.files.input2) "fastqc/${tnpair.tumour.name}.2")
    (ln (check tnpair.normal.files.input1) "fastqc/${tnpair.normal.name}.1")
    (ln (check tnpair.normal.files.input2) "fastqc/${tnpair.normal.name}.2")
  ];

in testNaming
