{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-tabix" { bgz = _: true; } input);

stage {

  name = "samtools-tabix";
  buildInputs = with pkgs; [ htslib ];
  buildCommand = ''
    ln -s ${input} input.bgz
    tabix input.bgz
    cp input.bgz.tbi $out
  '';
}
