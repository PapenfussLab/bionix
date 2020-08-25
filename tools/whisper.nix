{ bionix }:

with bionix;

rec {
  index = callBionixE ./whisper-index.nix;
  align = callBionixE ./whisper-align.nix;
}
