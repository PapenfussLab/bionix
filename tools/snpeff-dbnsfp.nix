{bionix
,nixpkgs
,dbnsfp
,flags ? ""}:

input:

with nixpkgs;
with bionix.types;

assert (matchFiletype "snpeff-dbnsfp" { vcf = _: true; } input);

stdenv.mkDerivation {
  name = "snpeff-dbnsfp";
  buildCommand = ''
    ln -s ${dbnsfp.db} dbNSFP.txt.gz
    ln -s ${dbnsfp.index} dbNSFP.txt.gz.tbi
    snpeff dbnsfp -db dbNSFP.txt.gz ${input} > $out
  '';
  buildInputs = [ snpeff ];
  passthru.filetype = input.filetype;
}
