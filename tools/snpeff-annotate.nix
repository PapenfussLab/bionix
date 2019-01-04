{bionix
,nixpkgs
,db
,flags ? ""}:

input:

with nixpkgs;
with bionix.types;

assert (matchFiletype "snpeff-annotate" { vcf = _: true; } input);

stdenv.mkDerivation {
  name = "snpeff-annotate";
  buildCommand = ''
    ln -s ${db} ${db.name}
    snpeff -nodownload -dataDir $TMPDIR ${db.name} ${input} > $out
  '';
  buildInputs = [ snpeff ];
  passthru.filetype = input.filetype;
}
