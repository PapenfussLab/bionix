{ bionix
, flags ? null
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "snap-index" { fa = _: true; } ref);

stage {
  name = "snap-index";
  buildInputs = [ bionix.snap.app ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    mkdir $out
    snap-aligner index ref.fa $out -t$NIX_BUILD_CORES ${optionalString (flags != null) flags}
  '';
  passthru.multicore = true;
}
