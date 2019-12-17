{ bionix, flags ? "", faidxAttrs ? {}, indexAttrs ? {}}:

inputs:

with bionix;
with lib;
with types;

let
  getref = matchFiletype "delly-call" { bam = x: x.ref; };
  refs = map getref inputs;
  ref = head refs;

  renameAndIndex = f:
    stage {
      name = "rename";
      buildCommand = ''
        mkdir $out
        ln -s ${f} $out/sample.bam
        ln -s ${samtools.index indexAttrs f} $out/sample.bam.bai
      '';
    };

in

assert (length (unique refs) == 1);

stage {
  name = "delly-call";
  buildInputs = with pkgs; [ delly bcftools ];
  buildCommand = ''
    export OMP_NUM_THREADS=$NIX_BUILD_CORES
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fa.fai
    delly call \
      -o delly.bcf \
      -g ref.fa \
      ${concatMapStringsSep " " (i: "${renameAndIndex i}/sample.bam") inputs} \
      ${flags}
    bcftools view delly.bcf > $out
  '';
  passthru = {
    multicore = true;
    filetype = filetype.vcf { ref = ref; };
  };
}
