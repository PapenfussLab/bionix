{ bionix
, nixpkgs
, flags ? null
, seed ? 42
}:

ref:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "bowtie-index" { fa = _: true; } ref);

stdenv.mkDerivation {
  name = "bowtie-index";
  buildInputs = [ bowtie2 ];
  buildCommand = ''
    mkdir $out
    bowtie2-build --seed ${toString seed} --threads $NIX_BUILD_CORES ${optionalString (flags != null) flags} ${ref} $out/ref
  '';
}
