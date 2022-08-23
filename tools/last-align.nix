{ bionix
, ref
, flags ? null
, trainFlags ? null
, indexAttrs ? { }
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;

let
  fq = f: matchFiletype "last-input" { fq = _: f; } f;
in
stage {
  name = "last-align";
  buildInputs = with pkgs; [ last ];
  buildCommand = ''
      last-train -P $NIX_BUILD_CORES \
        ${optionalString (trainFlags != null) flags} \
        ${bionix.lastal.index indexAttrs ref}/index \
        ${fq input1} ${optionalString (input2 != null) (fq input2)}\
        > train
      lastal -P $NIX_BUILD_CORES \
        ${optionalString (flags != null) flags} \
        ${bionix.lastal.index indexAttrs ref}/index \
        -p train \
        ${fq input1} ${optionalString (input2 != null) (fq input2)}\
        > tmp
    cp tmp $out
 '';
  passthru.multicore = true;
}
