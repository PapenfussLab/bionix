{bionix
,dbnsfp
,flags ? ""}:

input:

with bionix;
with types;

assert (matchFiletype "snpeff-dbnsfp" { vcf = _: true; } input);

stage {
  name = "snpeff-dbnsfp";
  buildCommand = ''
    ln -s ${dbnsfp.db} dbNSFP.txt.gz
    ln -s ${dbnsfp.index} dbNSFP.txt.gz.tbi
    snpsift dbnsfp -db dbNSFP.txt.gz ${input} > $out
  '';
  buildInputs = [ snpeff ];
  passthru.filetype = input.filetype;
}
