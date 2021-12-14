{ bionix
, flags ? null
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "subread-index" { fa = _: true; } ref);

stage {
  name = "subread-index";
  buildInputs = with pkgs; [ subread ];
  buildCommand = ''
    mkdir $out
    subread-buildindex ${optionalString (flags != null) flags} -o $out/ref ${ref}
  '';
}
