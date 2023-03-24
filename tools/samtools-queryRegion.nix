{ bionix
, flags ? null
, regions ? []
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "samtools-queryRegion" { fa = _: true; } input);

stage {

  name = "samtools-faidx";
  buildInputs = with pkgs; [ samtools ];
  buildCommand = ''
    ln -s ${input} input.fasta
    samtools faidx ${optionalString (flags != null) flags} input.fasta ${concatStringsSep " " regions} > $out
  '';
}
