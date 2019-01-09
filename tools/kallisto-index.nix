{bionix
, kmerSize ? 31
, unique ? false}:

with bionix;
with lib;
with types;

assert (kmerSize > 1);

input:

assert (matchFiletype input { fa = _: true; } input);

stage {
  name = "kallisto-index";
  buildInputs = with pkgs; [ kallisto ];
  buildCommand = ''
    kallisto index -k ${toString kmerSize} ${optionalString unique "--make-unique"} -i $out ${input}
  '';
}
