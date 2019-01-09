{ bionix
, flags ? null
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-faidx" { fa = _: true; } input);

stage {

  name = "samtools-faidx";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    ln -s ${input} input.fasta
    samtools faidx ${optionalString (flags != null) flags} input.fasta
    cp input.fasta.fai $out
  '';
}
