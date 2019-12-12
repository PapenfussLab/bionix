{ bionix
, faidxAttrs ? {}
, indexAttrs ? {}
, flags ? ""}:

with bionix;
with lib;
with types;

inputs:

let
  getref = f: matchFiletype "octopus-call" { bam = {ref, ...}: ref; cram = {ref, ...}: ref;} f;
  refs = map getref inputs;
  ref = head refs;

in

assert (length (unique refs) == 1);

stage {
  name = "octopus-call";
  buildInputs = with pkgs; [ octopus-caller ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fai
    ${concatMapStringsSep "\n" (i: ''
      ln -s ${i} $(basename ${i}).bam
      ln -s ${samtools.index indexAttrs i} $(basename ${i}).bai
    '') inputs}
    octopus -R ref.fa -I *.bam -o $out \
      --threads=$NIX_BUILD_CORES \
      ${flags}
  '';
  passthru.filtype = filetype.vcf {ref = ref;};
  passthru.multicore = true;
}
