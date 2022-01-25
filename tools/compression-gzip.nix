{ bionix, block ? false, flags ? "" }:

with bionix;
with types;

input:

stage {
  name = "gzip";
  buildInputs = with pkgs; [ (if block then htslib else gzip) ];
  buildCommand =
    if block then
      "bgzip -c ${flags} ${input} > $out"
    else
      "gzip ${flags} < ${input} > $out";
  passthru.filetype = (if block then filetype.bgz else filetype.gz) input.filetype;
}
