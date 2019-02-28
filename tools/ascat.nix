{bionix}:

with bionix;

{
  app = pkgs.callPackage ./ascat-app.nix {};
  gccorrect = callBionixE ./ascat-gccorrect.nix;
  callCNV = callBionixE ./ascat-callCNV.nix;
}
