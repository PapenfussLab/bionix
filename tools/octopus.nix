{ bionix }:

with bionix;

{
  /* Call variants for a population
    Type: { fast :: bool, very-fast :: bool, max-genotypes :: int, targets :: FilePath + [string], ... } -> [bam] -> vcf
  */
  call = callBionixE ./octopus-call.nix;

  /* Call somatic variants
    Type: { fast :: bool, very-fast :: bool, max-genotypes :: int, targets :: FilePath + [string], ... } -> { normal :: bam, tumours :: [bam] } -> vcf
  */
  callSomatic = callBionixE ./octopus-callSomatic.nix;
}
