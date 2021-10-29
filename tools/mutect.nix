{ bionix }:

with bionix;

{
  app = pkgs.callPackage ./mutect-app.nix { };

  /* Call somatic variants with mutect
    Type: { cosmic, dbsnp, ... } -> { normal :: bam, tumour :: bam } -> vcf
  */
  call = callBionixE ./mutect-call.nix;
}
