{ bionix, bwaIndexAttrs ? { }, faidxAttrs ? { }, indexAttrs ? { }
, assemblyAttrs ? { }, collectMetricsAttrs ? { }
, softClipsToSplitReadsAttrs ? { flags = "REALIGN_ENTIRE_READ=true"; }
, config ? null, heapSize ? "4g", shards ? 10 }:

with bionix;
with lib;
with types;
with gridss;

inputs:

let
  getref = matchFiletype "gridss-variants" { bam = x: x.ref; };
  ref = getref (head inputs);
  sorted = matchFileSorting "gridss-variants" { coord = _: true; };
  homoRef = length (unique (map getref inputs)) == 1;

  linkInput = f: attrs: input: ''
    BASENAME=$(basename ${input})
    WRKDIR="''${BASENAME}.gridss.working"
    if [[ ! -e $WRKDIR ]] ; then
      mkdir $WRKDIR
    fi
    for f in ${f attrs input}/* ; do
      ln -s $f $WRKDIR/$BASENAME.''${f##*.}
    done
  '';

  linkSV = input: ''
    BASENAME=$(basename ${input})
    WRKDIR="''${BASENAME}.gridss.working"
    if [[ ! -e $WRKDIR ]] ; then
      mkdir $WRKDIR
    fi
    ln -s ${input} $WRKDIR/$BASENAME.sv.bam
    ln -s ${bionix.samtools.index indexAttrs input} $WRKDIR/$BASENAME.sv.bai
  '';

  assembly = bionix.samtools.sort { }
    (softClipsToSplitReads softClipsToSplitReadsAttrs
      (bionix.gridss.shardedAssemble shards assemblyAttrs inputs));

  mkLinks = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    ${concatMapStringsSep "\n" (linkSV) inputs}
    ${linkSV assembly}
    ${concatMapStringsSep "\n" (linkInput collectMetrics collectMetricsAttrs)
    inputs}
    ${linkInput collectMetrics collectMetricsAttrs assembly}
  '';

in assert (all sorted inputs);
assert (homoRef);

rec {
  identify = stage rec {
    name = "gridss-identifyVariants";
    buildInputs = with pkgs; [ jre samtools ];
    buildCommand = mkLinks + ''
      java -Xmx${heapSize} -Dsamjdk.create_index=true \
        -cp ${jar} gridss.IdentifyVariants \
        VERBOSITY=WARNING \
        REFERENCE_SEQUENCE=ref.fa \
        ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
        ASSEMBLY=${assembly} \
        WORKER_THREADS=$NIX_BUILD_CORES \
        OUTPUT_VCF=out.vcf \
        ${
          optionalString (config != null)
          ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)
        } \
        WORKING_DIR=$TMPDIR/ \
        TMP_DIR=$TMPDIR/

      mv out.vcf $out
    '';
    passthru = {
      filetype = filetype.vcf { ref = ref; };
      gridss.assembly = assembly;
      multicore = true;
    };
  };

  annotate = stage rec {
    name = "gridss-annotateVariants";
    buildInputs = with pkgs; [ jre ];
    buildCommand = mkLinks + ''
      ln -s ${
        bionix.gridss.identifyVariants {
          inherit bwaIndexAttrs faidxAttrs indexAttrs assemblyAttrs
            collectMetricsAttrs softClipsToSplitReadsAttrs config;
        } inputs
      } input.vcf
      java -Xmx${heapSize} -Dsamjdk.create_index=true \
        -cp ${jar} gridss.AnnotateVariants \
        VERBOSITY=WARNING \
        REFERENCE_SEQUENCE=ref.fa \
        ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
        ASSEMBLY=${assembly} \
        INPUT_VCF=input.vcf \
        OUTPUT_VCF=out.vcf \
        WORKER_THREADS=$NIX_BUILD_CORES \
        WORKING_DIR=$TMPDIR/ \
        ${
          optionalString (config != null)
          ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)
        } \
        TMP_DIR=$TMPDIR/

      mv out.vcf $out
    '';
    passthru = {
      filetype = filetype.vcf { ref = ref; };
      gridss.assembly = assembly;
      multicore = true;
    };
  };

  annotateAndAssemble = stage rec {
    name = "gridss-annotateVariants-assemble";
    buildCommand = ''
      mkdir $out
      ln -s ${annotate} $out/gridss.vcf
      ln -s ${assembly} $out/gridss.bam
      ln -s ${bionix.samtools.index indexAttrs assembly} $out/gridss.bam.bai
    '';
  };
}
