{ bionix
, flags ? null
}:

gtf:

with bionix;
with lib;
with types;

stage {
  name = "hisat2-extractExons";
  buildInputs = with pkgs; [ hisat2 ];
  buildCommand = ''
    hisat2_extract_exons.py ${gtf} > $out
  '';
}
