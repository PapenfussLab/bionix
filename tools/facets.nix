{bionix}:

with bionix;

{
  app = lib.callPackageWith (pkgs // pkgs.rPackages) ./facets-app.nix {};
  callCNV = callBionix ./facets-call.nix;
}
