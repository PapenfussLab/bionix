{ bionix
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, indexAttrs ? {}
, collectMetricsAttrs ? {}
, flags ? null
, config ? null
, heapSize ? "31g"
, workdirs ? []
, jobIndex ? null
, jobNodes ? null
}:

with bionix;
with lib;
with types;

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
    ln -s ${input} $WRKDIR/$BASENAME.sv.bam
    ln -s  ${bionix.samtools.index indexAttrs input} $WRKDIR/$BASENAME.sv.bai
    for f in ${bionix.gridss.collectMetrics collectMetricsAttrs input}/* ; do
      ln -s $f $WRKDIR/$BASENAME.''${f##*.}
    done
  '';
in

assert (all sorted inputs);
assert (homoRef);

stage rec {
  name = "gridss-assemble";
  buildInputs = with pkgs; [ jre bwa rsync ];
  outputs = [ "out" "work" ];
  buildCommand = ''
    TMPDIR=$(pwd)
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    ${concatMapStringsSep "\n" linkInput inputs}
    ${concatMapStringsSep "\n" (w: "rsync -a --ignore-existing ${w}/ ./") workdirs}
    java -Xmx${heapSize} -Dsamjdk.create_index=true \
      -cp ${bionix.gridss.jar} gridss.AssembleBreakends \
      VERBOSITY=WARNING \
      REFERENCE_SEQUENCE=ref.fa \
      ${concatMapStringsSep " " (i: "INPUT='${i}'") inputs} \
      WORKER_THREADS=$NIX_BUILD_CORES \
      OUTPUT=$out \
      ${optionalString (config != null) ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)} \
      WORKING_DIR=$TMPDIR/ \
      TMP_DIR=$TMPDIR/ \
      ${optionalString (jobIndex != null) "JOB_INDEX=${toString jobIndex}"} \
      ${optionalString (jobIndex != null) "JOB_NODES=${toString jobNodes}"} \
      ${optionalString (flags != null) flags}
    rm -rf tmp
    touch $out
    cp -r $TMPDIR $work
    chmod u+rwX -R $work
  '';
  passthru.filetype = filetype.bam { ref = ref; sorting = sort.name {}; };
  passthru.multicore = true;
}
