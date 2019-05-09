{bionix}:

with bionix;

{
  app = pkgs.callPackage ./battenberg-app.nix {};
  callCNV = callBionixE ./battenberg-call.nix;
}
