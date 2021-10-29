{ bionix
, gtf ? null
, flags ? null
, extractSpliceSitesAttrs ? { }
, extractExonsAttrs ? { }
}:

ref:

with bionix;
with lib;
with types;

assert (matchFiletype "hisat2-index" { fa = _: true; } ref);

stage {
  name = "hisat2-index";
  buildInputs = with pkgs; [ hisat2 ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    mkdir $out
    hisat2-build -p $NIX_BUILD_CORES ${optionalString (flags != null) flags} \
      ${optionalString (gtf != null) "--ss ${hisat2.extractSpliceSites extractSpliceSitesAttrs gtf} --exon ${hisat2.extractExons extractExonsAttrs gtf}"} \
      ref.fa $out/index
  '';
  passthru.multicore = true;
}
