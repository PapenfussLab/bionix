{ bionix
, ref
, gc
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

{tumour, normal, gender}:

with bionix;
with lib;
with types;

stage rec {
  name = "ascat-callCNV";
  buildInputs = with pkgs; [ ascat.app ];
  buildCommand = ''
    mkdir $out
    ln -s ${tumour} tumour.bam
    ln -s ${bionix.samtools.index bamIndexAttrs tumour} tumour.bai
    ln -s ${normal} normal.bam
    ln -s ${bionix.samtools.index bamIndexAttrs normal} normal.bai
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx indexAttrs ref} ref.fa.fai
    ascat.pl \
      -outdir $out \
      -tumour tumour.bam \
      -normal normal.bam \
      -reference ref.fa \
      -snp_gc ${gc} \
      -gender ${gender} \
      -genderChr Y \
      -protocol WGS \
      -cpus $NIX_BUILD_CORES
  '';
  passthru.multicore = true;
}
