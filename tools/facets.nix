{ bionix }:

with bionix;

{
  app = lib.callPackageWith (pkgs // pkgs.rPackages) ./facets-app.nix { };

  /* Call CNVs
    Type: callCnv :: {...} -> {vcf, bams :: [bams]} -> CNVs
  */
  callCNV = callBionixE ./facets-call.nix;
}
