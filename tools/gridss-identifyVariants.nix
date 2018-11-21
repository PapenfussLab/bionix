{ bionix
, nixpkgs
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, indexAttrs ? {}
, assemblyAttrs ? {}
, collectMetricsAttrs ? {}
, softClipsToSplitReadsAttrs ? {}
, flags ? null
, config ? null
}:

with nixpkgs;
with lib;
with bionix.types;
with bionix.gridss;

inputs:

let
  getref = matchFiletype "gridss-identifyVariants" { bam = x: x.ref; };
  ref = getref (head inputs);
  sorted = matchFileSorting "gridss-identifyVariants" { coord = _: true; };
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

  linkSV = input: ''
    BASENAME=$(basename ${input})
    WRKDIR="''${BASENAME}.gridss.working"
    if [[ ! -e $WRKDIR ]] ; then
      mkdir $WRKDIR
    fi
    ln -s ${input} $WRKDIR/$BASENAME.sv.bam
    ln -s ${bionix.samtools.index indexAttrs input} $WRKDIR/$BASENAME.sv.bai
  '';

  assembly = bionix.samtools.sort {} (softClipsToSplitReads softClipsToSplitReadsAttrs (bionix.samtools.sort { nameSort = true;} (bionix.gridss.assemble assemblyAttrs inputs)));
in

assert (all sorted inputs);
assert (homoRef);

stdenv.mkDerivation rec {
  name = "gridss-identifyVariants";
  buildInputs = [ jre samtools ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    ${concatMapStringsSep "\n" (linkSV) inputs}
    ${linkSV assembly}
    ${concatMapStringsSep "\n" (linkInput collectMetrics collectMetricsAttrs) inputs}
    ${linkInput collectMetrics collectMetricsAttrs assembly}
	  java -Xmx4g -Dsamjdk.create_index=true \
      -cp ${jar} gridss.IdentifyVariants \
      REFERENCE_SEQUENCE=ref.fa \
      ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
      ASSEMBLY=${assembly} \
      OUTPUT_VCF=out.vcf \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + gridssConfig config)} \
      WORKING_DIR=$TMPDIR/ \
      TMP_DIR=$TMPDIR/

    mv out.vcf $out
    '';
  passthru = {
    filetype = filetype.vcf { ref = ref; };
    gridss.assembly = assembly;
  };
}
