{ bionix
, dbnsfp
, heapSize ? "31g"
, fields ? [ ]
, flags ? ""
}:

input:

with bionix;
with lib;
with types;

assert (matchFiletype "snpeff-dbnsfp" { vcf = _: true; } input);

stage {
  name = "snpeff-dbnsfp";
  JAVA_TOOL_OPTIONS = "-Xmx${heapSize}";
  buildCommand = ''
    ln -s ${dbnsfp.db} dbNSFP.txt.gz
    ln -s ${dbnsfp.index} dbNSFP.txt.gz.tbi
    snpsift dbnsfp -db dbNSFP.txt.gz ${optionalString (length fields > 0) "-f ${concatStringsSep "," fields}"} ${flags} ${input} > $out
  '';
  buildInputs = with pkgs; [ snpeff ];
  passthru.filetype = input.filetype;
}
