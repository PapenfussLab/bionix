{ bionix
, flags ? null
, altRegex ? "^>.*_alt$"
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "bwa-index" { fa = _: true; } ref);

stage {
  name = "bwa-index2";
  buildInputs = with pkgs; [ bionix.bwa.app2 ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    bwa-mem2 index ${optionalString (flags != null) flags} ref.fa
    mkdir $out
    mv ref.fa.* $out
    grep -P '${altRegex}' ref.fa | tr -d '^>' > $out/idxbase.alt || true
  '';
  stripStorePaths = false;
}
