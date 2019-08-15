{ bionix }:

with bionix;
with pkgs;

{
  app = callPackage ./snver-app.nix {};

  call = callBionix ./snver-call.nix;
}
