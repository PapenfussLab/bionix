{ bionix
, nixpkgs
, flags ? null
}:

ref:

with nixpkgs;
with lib;

stdenv.mkDerivation {
  name = "bwa-index";
  buildInputs = [ bwa ];
  buildCommand = ''
    ln -s ${ref.seq} ref.fa
    bwa index ${optionalString (flags != null) flags} ref.fa
    mkdir $out
    mv ref.fa.* $out
    grep '^>[^ \t]*_alt$' ref.fa | tr -d '^>' > $out/idxbase.alt || true
  '';
}
