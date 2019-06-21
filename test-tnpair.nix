with import <bionix> {};
with lib;

let
  fetchlocal = path: pkgs.stdenv.mkDerivation {
    name = baseNameOf path;
    buildCommand = "ln -s ${path} $out";
  };
  fetchfq = attrs: types.tagFiletype (types.filetype.fq {}) (fetchlocal attrs);
  fetchfa = attrs: types.tagFiletype (types.filetype.fa {}) (fetchlocal attrs);

  ref = fetchfa ./examples/ref.fa;

  alignWithRG = rg: bwa.align { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};
  sort = samtools.sort {};
  nameSort = samtools.sort {nameSort = true;};
  flagstat = samtools.flagstat {};
  check-fastqc = fastqc.check {};
  check-fastp = fastp.check {};
  callVariants = strelka.callSomatic {};
  markdup = samtools.markdup {};
  fixmate = samtools.fixmate {};

  tnpair = {
    tumour = {name = "mysample1"; files = {
        input1 = fetchfq ./examples/sample1-1.fq;
        input2 = fetchfq ./examples/sample1-2.fq;
      };
    };
    normal = {name = "mysample2"; files = {
        input1 = fetchfq ./examples/sample2-1.fq;
        input2 = fetchfq ./examples/sample2-2.fq;
      };
    };
  };

  processPair = { tumour, normal }: rec {
    alignments = mapAttrs (_: x: markdup (sort (fixmate (alignWithRG x.name x.files)))) { inherit normal tumour; };
    variants = callVariants alignments;
    glvariants = strelka.call {} (builtins.attrValues alignments);
    platypusVars = platypus.call {} (builtins.attrValues alignments);
    shards = pipe [
      (shard.fastQPair 2)
      (map (bwa.align {inherit ref;}))
    ] normal.files;
  };

  tnpairResult = processPair tnpair;

  cnvkitResults = rec {
    cnvs = cnvkit.callCNV {} (with tnpairResult.alignments; { normals = [ normal ]; tumours = [ tumour ];});
    plot = cnvkit.scatterPlot {} cnvs;
  };

  testNaming = linkDrv [
    (ln (facets.callCNV {} {vcf = tnpairResult.platypusVars; bams = with tnpairResult.alignments; [ normal tumour ];}) "facets")
    (ln cnvkitResults.cnvs "cnvkit")
    (ln cnvkitResults.plot "cnvkit.pdf")
    (ln tnpairResult.variants "strelka")
    (ln tnpairResult.glvariants "strelka-gl")
    (ln tnpairResult.variants.indels "strelka.indels.vcf")
    (ln tnpairResult.variants.snvs "strelka.snvs.vcf")
    (ln tnpairResult.glvariants.variants "strelka.gl.vcf")
    (ln (bowtie.align {inherit ref;} tnpair.normal.files) "alignments/bowtie-normal.bam")
    (ln (bwa.mem {inherit ref;} tnpair.normal.files) "alignments/bwa-mem.bam")
    (ln (bwa.mem2 {inherit ref;} tnpair.normal.files) "alignments/bwa-mem2.bam")
    (ln (minimap2.align {inherit ref; preset = "sr"; } tnpair.normal.files) "alignments/minimap2-normal.bam")
    (ln (snap.align {inherit ref; } tnpair.normal.files) "alignments/snap-normal.bam")
    (ln (gridss.callVariants {} (with tnpairResult.alignments; [normal tumour])) "gridss")
    (ln (gridss.call (with tnpairResult.alignments; [normal tumour])) "gridss2")
    (ln (gridss.callAndAssemble (with tnpairResult.alignments; [normal tumour])) "gridss3")
    (ln (samtools.merge {} tnpairResult.shards) "alignments/merged-shards.bam")
    (ln (samtools.merge {} [tnpairResult.alignments.tumour tnpairResult.alignments.normal]) "alignments/merged.bam")
    (ln (samtools.merge {} [(nameSort tnpairResult.alignments.tumour) (nameSort tnpairResult.alignments.normal)]) "alignments/merged-namesorted.bam")
    (ln (samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.tumour)) "alignments/${tnpair.tumour.name}.cram")
    #(ln (samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.normal)) "alignments/${tnpair.normal.name}.cram")
    (ln (flagstat tnpairResult.alignments.tumour) "alignments/${tnpair.tumour.name}.flagstat")
    #(ln (flagstat tnpairResult.alignments.normal) "alignments/${tnpair.normal.name}.flagstat")
    (ln (check-fastqc tnpair.tumour.files.input1) "fastqc/${tnpair.tumour.name}.1")
    #(ln (check-fastqc tnpair.normal.files.input1) "fastqc/${tnpair.normal.name}.1")
    #(ln (check-fastqc tnpair.normal.files.input2) "fastqc/${tnpair.normal.name}.2")
    #(ln (check-fastqc tnpair.tumour.files.input2) "fastqc/${tnpair.tumour.name}.2")
    (ln (check-fastp tnpair.tumour.files) "fastp/${tnpair.tumour.name}")
  ];

in testNaming
