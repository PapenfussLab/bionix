{ bionix
, ref
, bamOutput ? true
, flags ? null
, preset
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;
with compression;

let
  fa = f: matchFiletype "minimap2-ref" { fa = _: f; } f;
  fq = f: matchFiletype "minimap2-input" { fq = _: f; gz = matchFiletype' "minimap2-input" { fq = _: f; }; } f;

in
stage {
  name = "minimap2-align";
  buildInputs = with pkgs; [ minimap2 bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    ln -s ${fa ref} ref.fa
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      cores=1
    fi
    minimap2 ${optionalString (flags != null) flags} -t $cores -ax ${preset} ref.fa ${fq input1} \
      ${optionalString (input2 != null) (fq input2)} \
      ${optionalString bamOutput "| samtools view -b"} \
      > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam { inherit ref; sorting = sort.none { }; } else filetype.sam { inherit ref; sorting = sort.none { }; };
  passthru.multicore = true;
}
