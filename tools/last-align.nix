{ bionix
, ref
, flags ? null
, indexAttrs ? {}
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;

let
  fa = f: matchFiletype "last-ref" { fa = _: f; } f;
  fq = f: matchFiletype "last-input" { fq = _: f; } f;
  fqfa = f: matchFiletype "last-input" { fq = _: f; fa = _: f; } f;
  Q = f: matchFiletype "last-input" { fq = _: 1; fa = _: 0; } f;
  parallel = f: matchFiletype "last-input" { fq = _: "parallel-fastq"; fa = _: "parallel-fasta"; } f;

in stage {
  name = "last-align";
  buildInputs = with pkgs; [ last ];
  buildCommand = ''
    ${if (input2 != null) then "fastq-interleave ${fq input1} ${fq input2}" else "cat ${fa input1}"} | \
      ${parallel input1} -j $NIX_BUILD_CORES \
      lastal -Q${toString (Q input1)} \
        ${optionalString (input2 != null) "-i1"} \
        ${optionalString (flags != null) flags} \
        ${bionix.lastal.index indexAttrs ref}/index \
        > tmp
  '' + (if input2 == null then ''
    cp tmp $out
  '' else ''
    last-pair-probs tmp > $out
  '');
  passthru.multicore = true;
}
