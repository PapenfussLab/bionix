{ bionix, flags ? "", faidxAttrs ? { }, indexAttrs ? { } }:

inputs:

with bionix;
with lib;
with types;

let
  getref = matchFiletype "lumpy-call" { bam = x: x.ref; };
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

in assert (length (unique refs) == 1);

stage {
  name = "lumpy";
  buildInputs = with pkgs; [ lumpy ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fa.fai
    lumpyexpress \
      ${
        concatMapStringsSep " " (i: "-B ${renameAndIndex i}/sample.bam") inputs
      } \
      ${flags} \
      -o $out
  '';
  passthru.filetype = filetype.vcf { ref = ref; };
}
