{ bionix
, flags ? null
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "bwa-index" { fa = _: true; } ref);

stage {
  name = "bwa-index";
  buildInputs = with pkgs; [ bwa ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    bwa index ${optionalString (flags != null) flags} ref.fa
    mkdir $out
    mv ref.fa.* $out
    grep '^>[^ \t]*_alt$' ref.fa | tr -d '^>' > $out/idxbase.alt || true
  '';
}
