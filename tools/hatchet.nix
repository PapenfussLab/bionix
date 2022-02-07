{ bionix }:

with bionix;

{
  app = pkgs.callPackage ./hatchet-app.nix { };
  call = callBionixE ./hatchet-call.nix;
}
