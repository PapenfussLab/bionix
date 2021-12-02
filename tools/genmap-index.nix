{ bionix }:

with bionix;
with lib.types;

let
  fa = f: matchFiletype "genmap-index" { fa = _: f; } f;
in

ref: stage {
  name = "genmap-index";
  buildInputs = with pkgs; [ genmap ];
  buildCommand = ''
    genmap index -F ${fa ref} -I $out
  '';
}
