{ bionix
, ref
, bamOutput ? true
, flags ? null
, indexAttrs ? {}
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;
with compression;

let
  fa = f: matchFiletype "star-ref" { fa = _: f; } f;
  fq = f: matchFiletype "star-input" { fq = _: f; gz = matchFiletype' "star-input" { fq = _: "<(gunzip < ${f})"; }; } f;

in stage {
  name = "star-align";
  buildInputs = with pkgs; [ star bc samtools ];
  buildCommand = ''
    ln -s ${fa ref} ref.fa
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      cores=1
    fi
    STAR ${optionalString (flags != null) flags} \
     --runThreadN $cores \
     --genomeDir ${star.index indexAttrs ref} \
     --readFilesIn ${fq input1} ${optionalString (input2 != null) (fq input2)}
     ${if bamOutput then "samtools view -b Aligned.out.sam > $out" else "cp Aligned.out.sam $out"}
  '';
  passthru.filetype = if bamOutput then filetype.bam {ref = ref; sorting = sort.none {};} else filetype.sam {ref = ref; sorting = sort.name {};};
  passthru.multicore = true;
}
