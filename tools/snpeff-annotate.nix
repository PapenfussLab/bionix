{ bionix
, db
, heapSize ? "31g"
, flags ? ""
}:

input:

with bionix;
with types;

assert (matchFiletype "snpeff-annotate" { vcf = _: true; } input);

stage {
  name = "snpeff-annotate";
  buildInputs = with pkgs; [ snpeff ];
  JAVA_TOOL_OPTIONS = "-Xmx${heapSize}";
  buildCommand = ''
    ln -s ${db} ${db.name}
    snpeff -nodownload -dataDir $TMPDIR ${db.name} ${input} > $out
  '';
  passthru.filetype = input.filetype;
}
