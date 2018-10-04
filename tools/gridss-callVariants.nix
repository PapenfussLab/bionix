{ bionix
, nixpkgs
, ref
, blacklist ? null
, bwaIndexAttrs ? {}
, faidxAttrs ? {}
, flags ? null
}:

with lib;

inputs:

stdenv.mkDerivation rec {
  name = "gridss-callVariants";
  buildInputs = [ jre R bwa ];
  jar = fetchurl {
    url = "https://github.com/PapenfussLab/gridss/releases/download/v2.0.0/gridss-2.0.0-gridss-jar-with-dependencies.jar";
    sha256 = "01srl3qvv060whqg1y1fpxjc5cwga5wscs1bmf1v3z87dignra7k";
  };
  buildCommand = ''
    ln -s ${ref.seq} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    for f in ${bionix.bwa.index bwaIndexAttrs ref}/*; do
      ln -s $f
    done
    mkdir $out
    java -ea -Xmx31g \
	    -Dreference_fasta="ref.fa" \
	    -Dsamjdk.create_index=true \
	    -Dsamjdk.use_async_io_read_samtools=true \
	    -Dsamjdk.use_async_io_write_samtools=true \
	    -Dsamjdk.use_async_io_write_tribble=true \
	    -Dgridss.gridss.output_to_temp_file=true \
	    -cp ${jar} gridss.CallVariants \
      WORKER_THREADS=$NIX_BUILD_CORES \
	    TMP_DIR=. \
	    WORKING_DIR=. \
	    REFERENCE_SEQUENCE="ref.fa" \
      ${concatMapStringsSep " " (i: "INPUT=\"${i}\"") inputs} \
	    OUTPUT="$out/gridss.vcf" \
	    ASSEMBLY="$out/gridss.bam" \
      ${optionalString (blacklist != null) ("BLACKLIST=" + blacklist)} \
      ${optionalString (flags != null) flags}
  '';
}
