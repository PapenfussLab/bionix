{ bionix }:

with bionix;

rec {
  app = pkgs.callPackage ./whisper-app.nix {};
  index = callBionixE ./whisper-index.nix;
  align = callBionixE ./whisper-align.nix;
}
