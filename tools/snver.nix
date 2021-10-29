{ bionix }:

with bionix;
with pkgs;

{
  app = callPackage ./snver-app.nix { };

  /* Call variants
    Type: { ploidy, ... } -> [bam] -> vcf
  */
  call = callBionix ./snver-call.nix;
}
