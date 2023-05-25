{bionix}:
with bionix;

{
  app = pkgs.callPackage ./aa-app.nix { };
  call = callBionixE ./aa-call.nix { };
}
