{ bionix
, ref
, bamOutput ? true
, flags ? null
, indexAttrs ? { }
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;
with compression;

let
  fa = f: matchFiletype "whisper-ref" { fa = _: f; } f;
  fq = f: matchFiletype "whisper-input" { fq = _: f; gz = matchFiletype' "whisper-input" { fq = _: f; }; } f;

in
stage {
  name = "whisper-mem";
  buildInputs = with pkgs; [ whisper ];
  buildCommand = ''
    ln -s ${fa ref} ref.fa
    for f in ${bionix.whisper.index indexAttrs ref}/* ; do
      ln -s $f
    done
    whisper ${optionalString (flags != null) flags} -t $NIX_BUILD_CORES \
      ${optionalString bamOutput "-store-BAM"} \
      ${optionalString (input2 != null) "-rp"} \
      -stdout \
      index ${fq input1} \
      ${optionalString (input2 != null) (fq input2)} \
      > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam { inherit ref; sorting = sort.none { }; } else filetype.sam { inherit ref; sorting = sort.name { }; };
  passthru.multicore = true;
  stripStorePaths = false;
}
