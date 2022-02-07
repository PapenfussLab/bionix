{ bionix }:

with bionix;

{
  app = pkgs.callPackage ./pizzly-app.nix { };
  call = callBionixE ./pizzly-call.nix;
}
