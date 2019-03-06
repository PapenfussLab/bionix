{bionix}:

with bionix;

{
  app = lib.callPackageWith (pkgs // pkgs.rPackages) ./facets-app.nix {};

  /* Call CNVs
  Type: callCnv :: {...} -> {vcf, bams :: [bams]} -> CNVs
  */
  callCNV = callBionix ./facets-call.nix;
}
