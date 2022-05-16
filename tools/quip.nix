{bionix}:
with bionix; {
  app = pkgs.callPackage ./quip-app.nix {};
  quip = callBionixE ./quip-quip.nix;
  unquip = callBionixE ./quip-unquip.nix;
}
