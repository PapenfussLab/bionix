{ bionix
, flags ? null
}:

gtf:

with bionix;
with lib;
with types;

stage {
  name = "hisat2-extractSpliceSites";
  buildInputs = with pkgs; [ hisat2 ];
  buildCommand = ''
    hisat2_extract_splice_sites.py ${gtf} > $out
  '';
}
