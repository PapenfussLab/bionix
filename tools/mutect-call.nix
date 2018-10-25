{bionix
, nixpkgs
, cosmic
, dbsnp}:

with nixpkgs;
with lib;

let
  inherit (bionix.types) matchFiletype;
  getVCFref = matchFiletype "mutect-call" {vcf = {ref}: ref;};
  getBAMref = matchFiletype "mutect-call" {bam = {ref, ...}: ref;};
  refs = map getVCFref [ cosmic dbsnp ];
  ref = head refs;
in

assert (length (unique refs) == 1);

{normal, tumour}:

assert (ref == getBAMref normal && ref == getBAMref tumour);

stdenv.mkDerivation {
  name = "mutect";
  buildInputs = [ bionix.mutect.app ];
  buildCommand = ''
    ln -s ${normal} normal.bam
    ln -s ${tumour} tumour.bam
    ln -s ${dbsnp} dbsnp.vcf
    ln -s ${cosmic} cosmic.vcf
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx {} ref} ref.fa.fai
    ln -s ${bionix.samtools.dict {} ref} ref.dict
    mutect --analysis_type MuTect \
      --reference_sequence ref.fa \
      --cosmic cosmic.vcf \
      --dbsnp dbsnp.vcf \
      --input_file:normal normal.bam \
      --input_file:tumour tumour.bam \
      --out $out
  '';
}
