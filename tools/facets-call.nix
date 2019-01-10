{bionix}:

{vcf, bams}:

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
    # Facets requires lexical sorting on the VCF files
    grep '^#' ${vcf} > input.vcf
    grep -v '^#' ${vcf} | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> input.vcf || true

    # Now actually run facets
    snp-pileup input.vcf $out ${concatStringsSep " " bams}
  '';
}
