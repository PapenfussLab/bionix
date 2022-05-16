{bionix}:
with bionix; {
  app = pkgs.callPackage ./quip-app.nix {};
  quip = callBionix ./quip-quip.nix;
  unquip = callBionix ./quip-unquip.nix;
}
