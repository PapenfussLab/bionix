{ bionix
, ref
, chrPrefix ? ""
, flags ? null
}:

snp:

with bionix;
with lib;
with types;

stage rec {
  name = "ascat-gccorrect";
  buildInputs = with pkgs; [ ascat.app gawk ];
  script = pkgs.writeText "convert.awk" ''
    BEGIN{
      FS = OFS = "\t"
    }
    /^#/{next}
    !loc[$1,$2]{
      print $3, "${chrPrefix}" $1, $2
      loc[$1,$2]++
    }
  '';
  buildCommand = ''
    awk -f ${script} ${snp} > snpPos.tsv
    mkdir splitPos splitGc splitGcLogs
    split --number=l/$NIX_BUILD_CORES -d snpPos.tsv splitPos/snpPos.
    ls -1 splitPos/ | xargs -n1 -P$NIX_BUILD_CORES -I '{}' sh -c 'ascatSnpPanelGcCorrections.pl ${ref} splitPos/{} > splitGc/{}'
    mv splitGc/snpPos.00 $out
    for f in splitGc/* ; do
      sed 1d $f >> $out
    done
  '';
  passthru.multicore = true;
}
