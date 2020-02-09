{ bionix
, flags ? null
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "last-index" { fa = _: true; } ref);

stage {
  name = "last-index";
  buildInputs = with pkgs; [ last ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    mkdir $out
    lastdb ${optionalString (flags != null) flags} $out/index ref.fa
  '';
}
