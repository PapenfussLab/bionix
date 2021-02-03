{ bionix }:

with bionix;

{
  /* Call structural variants
  Type: call { ... } -> [bam] -> vcf
  */
  call = callBionixE ./delly-call.nix;
}
