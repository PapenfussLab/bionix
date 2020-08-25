{ bionix
, flags ? null
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "whisper-index" { fa = _: true; } ref);

stage {
  name = "whisper-index";
  buildInputs = with pkgs; [ whisper ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    mkdir $out
    whisper-index ${optionalString (flags != null) flags} index ref.fa $out $TMPDIR
  '';
}
