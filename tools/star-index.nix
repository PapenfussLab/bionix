{ bionix
, gtf
, flags ? null
, extractSpliceSitesAttrs ? {}
, extractExonsAttrs ? {}
, overhang ? 100
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "star-index" { fa = _: true; } ref);

stage {
  name = "star-index";
  buildInputs = with pkgs; [ star ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    mkdir $out
    STAR --runMode genomeGenerate \
      --runThreadN $NIX_BUILD_CORES \
      --sjdbGTFfile ${gtf} \
      --genomeDir $out \
      --genomeFastaFiles ${ref} \
      --sjdbOverhang ${toString overhang} \
      ${optionalString (flags != null) flags}
  '';
  passthru.multicore = true;
  stripStorePaths = false;
}
