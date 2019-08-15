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
  fa = f: matchFiletype "hisat2-ref" { fa = _: f; } f;
  fq = f: matchFiletype "hisat2-input" { fq = _: f; gz = matchFiletype' "hisat2-input" { fq = _: "<(gunzip < ${f})"; }; } f;

in stage {
  name = "hisat2-align";
  buildInputs = with pkgs; [ hisat2 bc samtools ];
  buildCommand = ''
    ln -s ${fa ref} ref.fa
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      cores=1
    fi
    hisat2 ${optionalString (flags != null) flags} -p $cores -x ${hisat2.index indexAttrs ref}/ \
      ${if input2 != null then "-1 ${fq input1} -2 ${fq input2}" else "-U ${fq input1}"} \
      ${optionalString bamOutput "| samtools view -b"} \
      | samtools sort -n \
      > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam {ref = ref; sorting = sort.name {};} else filetype.sam {ref = ref; sorting = sort.name {};};
  passthru.multicore = true;
}
