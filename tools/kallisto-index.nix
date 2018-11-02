{bionix
, nixpkgs
, kmerSize ? 31
, unique ? false}:

with nixpkgs;
with lib;
with bionix.types;

assert (kmerSize > 1);

input:

assert (matchFiletype input { fa = _: true; } input);

stdenv.mkDerivation {
  name = "kallisto-index";
  buildInputs = [ kallisto ];
  buildCommand = ''
    kallisto index -k ${toString kmerSize} ${optionalString unique "--make-unique"} -i $out ${input}
  '';
}
