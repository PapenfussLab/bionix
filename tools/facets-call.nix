{ bionix, vcf }:

bams:

with bionix;
with types;
with lib;

assert (matchFiletype "facets-call-vcf" { vcf = _: true; } vcf);
assert (all (matchFiletype "facets-call-bam" { bam = _: true; }) bams);
assert (all (matchFileSorting "facets-call-bam" { coord = _: true; }) bams);

stage {
  name = "facets";
  buildInputs = [ facets.app ];
  buildCommand = ''
    snp-pileup ${vcf} $out ${concatStringsSep " " bams}
  '';
}
