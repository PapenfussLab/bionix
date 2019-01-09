{bionix
,db
,flags ? ""}:

input:

with bionix;
with types;

assert (matchFiletype "snpeff-annotate" { vcf = _: true; } input);

stage {
  name = "snpeff-annotate";
  buildInputs = with pkgs; [ snpeff ];
  buildCommand = ''
    ln -s ${db} ${db.name}
    snpeff -nodownload -dataDir $TMPDIR ${db.name} ${input} > $out
  '';
  passthru.filetype = input.filetype;
}
