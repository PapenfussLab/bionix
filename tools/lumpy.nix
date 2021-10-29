{ bionix }:

with bionix;

{
  /* Call structural variants
    Type: { ... } -> [bam] -> vcf
  */
  call = callBionixE ./lumpy-call.nix;
}
