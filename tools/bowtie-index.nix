{ bionix
, flags ? null
, seed ? 42
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "bowtie-index" { fa = _: true; } ref);

stage {
  name = "bowtie-index";
  buildInputs = with pkgs; [ bowtie2 ];
  buildCommand = ''
    mkdir $out
    bowtie2-build --seed ${toString seed} --threads $NIX_BUILD_CORES ${optionalString (flags != null) flags} ${ref} $out/ref
  '';
  passthru.multicore = true;
}
