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
    octopusSomatic = octopus.callSomatic {} {inherit (alignments) normal; tumours = [ alignments.tumour ];};
    glvariants = strelka.call {} (builtins.attrValues alignments);
    platypusVars = platypus.call {} (builtins.attrValues alignments);
    octopusVars = octopus.call {} (builtins.attrValues alignments);
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

  alignments = {
    "bowtie-normal.bam" = bowtie.align {inherit ref;} tnpair.normal.files;
    "bwa-mem.bam" = bwa.mem {inherit ref;} tnpair.normal.files;
    "bwa-mem2.bam" = bwa.mem2 {inherit ref;} tnpair.normal.files;
    "minimap2-normal.bam" = minimap2.align {inherit ref; preset = "sr"; } tnpair.normal.files;
    "snap-normal.bam" = snap.align {inherit ref; } tnpair.normal.files;
    "${tnpair.tumour.name}.flagstat" = flagstat tnpairResult.alignments.tumour;
  };

  testNaming = linkOutputs {
    facets = facets.callCNV {} {vcf = tnpairResult.platypusVars; bams = with tnpairResult.alignments; [ normal tumour ];};
    cnvkit = cnvkitResults.cnvs;
    "cnvkit.pdf" = cnvkitResults.plot;
    "octopus.vcf" = tnpairResult.octopusVars;
    "octopus-somatic.vcf" = tnpairResult.octopusSomatic;
    strelka = tnpairResult.variants;
    strelka-gl = tnpairResult.glvariants;
    strelka-indels = tnpairResult.variants.indels;
    "strelka.snvs.vcf" = tnpairResult.variants.snvs;
    "strelka.gl.vcf" = tnpairResult.glvariants.variants;
    gridss = gridss.callVariants {} (with tnpairResult.alignments; [normal tumour]);
    gridss2 = gridss.call (with tnpairResult.alignments; [normal tumour]);
    gridss3 = gridss.callAndAssemble (with tnpairResult.alignments; [normal tumour]);
    "merged-shards.bam" = samtools.merge {} tnpairResult.shards;
    "merged.bam" = samtools.merge {} [tnpairResult.alignments.tumour tnpairResult.alignments.normal];
    "merged-namesorted.bam" = samtools.merge {} [(nameSort tnpairResult.alignments.tumour) (nameSort tnpairResult.alignments.normal)];
    "${tnpair.tumour.name}.cram" = samtools.view { outfmt = types.toCram; } (tnpairResult.alignments.tumour);
    "${tnpair.tumour.name}.fastqc.1" = check-fastqc tnpair.tumour.files.input1;
    "${tnpair.tumour.name}.fastp" = check-fastp tnpair.tumour.files;
    inherit alignments;
  };

in testNaming
