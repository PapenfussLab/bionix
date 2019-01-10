{bionix}:

with bionix;

{
  app = pkgs.callPackage ./mutect-app.nix {};
  call = callBionixE ./mutect-call.nix;
}
