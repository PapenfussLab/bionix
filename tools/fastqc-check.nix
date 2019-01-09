{ bionix
, flags ? null
}:

with bionix;
with lib;

input:

stage {
  name = "fastqc-check";
  buildInputs = [ bionix.fastqc.fastqc ];
  buildCommand = ''
    mkdir $out
    fastqc \
      -o $out \
      ${optionalString (flags != null) flags} \
      ${input}
  '';
}
