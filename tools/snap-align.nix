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
  fa = f: matchFiletype "snap-ref" { fa = _: f; } f;
  fq = f: matchFiletype "snap-input" { fq = _: f; gz = matchFiletype' "snap-input" { fq = _: f; }; } f;

in
stage {
  name = "snap-align";
  buildInputs = with pkgs; [ bionix.snap.app bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    for f in /* ; do
      ln -s $f
    done
    snap-aligner ${if input2 == null then "single" else "paired"} ${bionix.snap.index indexAttrs ref} ${fq input1} \
      ${optionalString (input2 != null) (fq input2)} \
      -o ${if bamOutput then "-bam" else "-sam" } - \
      -t $NIX_BUILD_CORES \
      ${optionalString (flags != null) flags} \
      | samtools sort -n > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam { inherit ref; sorting = sort.none { }; } else filetype.sam { inherit ref; sorting = sort.name { }; };
  passthru.multicore = true;
  stripStorePaths = false;
}
