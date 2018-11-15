{ bionix
, nixpkgs
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, collectMetricsAttrs ? {}
, extractSVReadsAttrs ? {}
, flags ? null
}:

with nixpkgs;
with lib;
with bionix.types;

inputs:

let
  getref = matchFiletype "gridss-assemble" { bam = x: x.ref; };
  ref = getref (head inputs);
  sorted = matchFileSorting "gridss-assemble" { coord = _: true; };
  homoRef = length (unique (map getref inputs)) == 1;

  linkInput = input: ''
    BASENAME=$(basename ${input})
    WRKDIR="''${BASENAME}.gridss.working"
    mkdir $WRKDIR
    for f in ${bionix.gridss.extractSVReads extractSVReadsAttrs input}/* ; do
      ln -s $f $WRKDIR/$BASENAME.''${f#*.}
    done
    for f in ${bionix.gridss.collectMetrics collectMetricsAttrs input}/* ; do
      ln -s $f $WRKDIR/$BASENAME.''${f#*.}
    done
  '';
in

assert (all sorted inputs);
assert (homoRef);

stdenv.mkDerivation rec {
  name = "gridss-assemble";
  buildInputs = [ jre bwa ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    ${concatMapStringsSep "\n" linkInput inputs}
	  java -Xmx31g -Dsamjdk.create_index=true \
      -cp ${bionix.gridss.jar} gridss.AssembleBreakends \
      REFERENCE_SEQUENCE=ref.fa \
      ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
      WORKER_THREADS=$NIX_BUILD_CORES \
      OUTPUT=$out \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + bionix.gridss.gridssConfig config)} \
      WORKING_DIR=$TMPDIR/ \
      TMP_DIR=$TMPDIR/
  '';
  passthru.filetype = filetype.bam { ref = ref; sorting = sort.coord {}; };
}
