{ bionix
, nixpkgs
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, assemblyAttrs ? {}
, extractSVReadsAttrs ? {}
, collectMetricsAttrs ? {}
, softClipsToSplitReadsAttrs ? {}
, identifyVariantsAttrs ? {}
, flags ? null
, config ? null
}:

with nixpkgs;
with lib;
with bionix.types;
with bionix.gridss;

inputs:

let
  getref = matchFiletype "gridss-annotateVariants" { bam = x: x.ref; };
  ref = getref (head inputs);
  sorted = matchFileSorting "gridss-annotateVariants" { coord = _: true; };
  homoRef = length (unique (map getref inputs)) == 1;

  linkInput = f: attrs: input: ''
    BASENAME=$(basename ${input})
    WRKDIR="''${BASENAME}.gridss.working"
    if [[ ! -e $WRKDIR ]] ; then
      mkdir $WRKDIR
    fi
    for f in ${f attrs input}/* ; do
      ln -s $f $WRKDIR/$BASENAME.''${f#*.}
    done
  '';

  assembly = bionix.samtools.sort {} (softClipsToSplitReads softClipsToSplitReadsAttrs (bionix.samtools.sort { nameSort = true;} (bionix.gridss.assemble assemblyAttrs inputs)));
in

assert (all sorted inputs);
assert (homoRef);

stdenv.mkDerivation rec {
  name = "gridss-identifyVariants";
  buildInputs = [ jre ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    ${concatMapStringsSep "\n" (linkInput extractSVReads extractSVReadsAttrs) inputs}
    ${concatMapStringsSep "\n" (linkInput collectMetrics collectMetricsAttrs) inputs}
    ${linkInput collectMetrics collectMetricsAttrs assembly}
    ASSBASE=$(basename ${assembly})
    ln -s ${assembly} $ASSBASE.gridss.working/$ASSBASE.sv.bam
    ln -s ${bionix.samtools.index {} assembly} $ASSBASE.gridss.working/$ASSBASE.sv.bai
    ln -s ${identifyVariants identifyVariantsAttrs inputs} input.vcf
	  java -Xmx4g -Dsamjdk.create_index=true \
      -cp ${jar} gridss.AnnotateVariants \
      REFERENCE_SEQUENCE=ref.fa \
      ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
      ASSEMBLY=${assembly} \
      INPUT_VCF=input.vcf \
      OUTPUT_VCF=out.vcf \
      WORKING_DIR=$TMPDIR/ \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + gridssConfig config)} \
      TMP_DIR=$TMPDIR/

    mv out.vcf $out
    '';
  passthru = {
    filetype = filetype.vcf { ref = ref; };
    gridss.assembly = assembly;
  };
}
