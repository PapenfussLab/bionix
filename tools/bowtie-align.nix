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
  fa = f: matchFiletype "bowtie2-ref" { fa = _: f; } f;
  fq = f: matchFiletype "bowtie2-input" { fq = _: f; gz = matchFiletype' "bowtie2-input" { fq = _: f; }; } f;

in stage {
  name = "bowtie2-align";
  buildInputs = with pkgs; [ bowtie2 bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      cores=1
    fi
    bowtie2 -x ${bionix.bowtie.index indexAttrs ref}/ref ${optionalString (flags != null) flags} --threads $cores \
      ${if input2 != null then "-1 " + fq input1 + " -2 " + fq input2 else "-U " + fq input1} \
      ${optionalString bamOutput "| samtools view -b"} \
      > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam {ref = ref; sorting = sort.name {};} else filetype.sam {ref = ref; sorting = sort.name {};};
  passthru.multicore = true;
}
