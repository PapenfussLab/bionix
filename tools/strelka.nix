{ bionix }:

with bionix;
with types;

{
  callSomatic = callBionixE ./strelka-callSomatic.nix;
  call = callBionixE ./strelka-call.nix;
  variants = drv: stage {
    name = "strelka-call-variants";
    buildCommand = ''
      ln -s ${drv}/variants/variants.vcf.gz $out
    '';
    passthru.filetype = filetype.gz (filetype.vcf {ref=ref;});
  };
  indels = drv: stage {
    name = "strelka-callVariants-indels";
    buildCommand = "ln -s ${drv}/variants/somatic.indels.vcf.gz $out";
    passthru.filetype = filetype.gz (filetype.vcf {ref = ref;});
  };
  snvs = drv: stage {
    name = "strelka-callVariants-snvs";
    buildCommand = "ln -s ${drv}/variants/somatic.snvs.vcf.gz $out";
    passthru.filetype = filetype.gz (filetype.vcf {ref = ref;});
  };
}
