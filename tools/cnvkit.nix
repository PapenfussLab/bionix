{bionix}:

with bionix;

{
app = lib.callPackageWith (pkgs // pkgs.pythonPackages) ./cnvkit-app.nix {};
  callCNV = callBionixE ./cnvkit-batch.nix;
}
