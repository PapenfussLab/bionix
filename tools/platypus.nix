{ bionix }:

with bionix;

{
  /* Call variants 
    Type: { ... } -> [bam] -> vcf
  */
  call = callBionixE ./platypus-callVariants.nix;
}
