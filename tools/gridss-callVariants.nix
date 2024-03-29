{ bionix
, blacklist ? null
, bwaIndexAttrs ? { }
, faidxAttrs ? { }
, flags ? null
, config ? null
, heapSize ? "31g"
}:

with bionix;
with lib;
with types;

inputs:

let
  getref = matchFiletype "gridss-callVariants" { bam = x: x.ref; };
  refs = map getref inputs;
  ref = head refs;
in

assert (length (unique refs) == 1);

stage rec {
  name = "gridss-callVariants";
  buildInputs = with pkgs; [ jre R bwa ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    mkdir $out
    java -ea -Xmx${heapSize} \
      -Dreference_fasta="ref.fa" \
      -Dsamjdk.create_index=true \
      -Dsamjdk.use_async_io_read_samtools=true \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=true \
      -Dgridss.gridss.output_to_temp_file=true \
      -cp ${bionix.gridss.jar} gridss.CallVariants \
      VERBOSITY=WARNING \
      WORKER_THREADS=$NIX_BUILD_CORES \
      TMP_DIR=. \
      WORKING_DIR=. \
      ${optionalString (config != null) ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)} \
      REFERENCE_SEQUENCE="ref.fa" \
      ${concatMapStringsSep " " (i: "INPUT=\"${i}\"") inputs} \
      OUTPUT="$out/gridss.vcf" \
      ASSEMBLY="$out/gridss.bam" \
      ${optionalString (blacklist != null) ("BLACKLIST=" + blacklist)} \
      ${optionalString (flags != null) flags}

    # The VCF index is non-deterministic
    rm $out/gridss.vcf.idx
  '';
  passthru.multicore = true;
}
