{ bionix }:

with bionix;
with types;

{
  /* Calls somatic variants
  Type: callSomatic :: {...} -> {tumour, normal} -> somatic results
  */
  callSomatic = callBionixE ./strelka-callSomatic.nix;
  /* Calls variants
  Type: call :: {...} -> [input] -> results
  */
  call = callBionixE ./strelka-call.nix;
  /* Extract VCF file from results
  Type: variants :: results -> vcf
  */
  variants =
    # result of call
    drv: stage {
    name = "strelka-call-variants";
    buildCommand = ''
      ln -s ${drv}/variants/variants.vcf $out
    '';
    passthru.filetype = filetype.vcf {ref=ref;};
  };
  /* Extract indels from somatic results
  Type: indels :: somatic results -> vcf
  */
  indels =
    # result of callSomatic
    drv: stage {
    name = "strelka-callVariants-indels";
    buildCommand = "ln -s ${drv}/variants/somatic.indels.vcf $out";
    passthru.filetype = filetype.vcf {ref = ref;};
  };
  /* Extract SNVs from somatic results
  Type: snvs :: somatic results -> vcf
  */
  snvs =
    # result of callSomatic
    drv: stage {
    name = "strelka-callVariants-snvs";
    buildCommand = "ln -s ${drv}/variants/somatic.snvs.vcf $out";
    passthru.filetype = filetype.vcf {ref = ref;};
  };
}
