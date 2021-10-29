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
  fa = f: matchFiletype "bwa-ref" { fa = _: f; } f;
  fq = f: matchFiletype "bwa-input" { fq = _: f; gz = matchFiletype' "bwa-input" { fq = _: f; }; } f;

in
stage {
  name = "bwa-mem";
  buildInputs = with pkgs; [ bwa bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    ln -s ${fa ref} ref.fa
    for f in ${bionix.bwa.index indexAttrs ref}/* ; do
      ln -s $f
    done
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      cores=1
    fi
    bwa mem ${optionalString (flags != null) flags} -t $cores ref.fa ${fq input1} \
      ${optionalString (input2 != null) (fq input2)} \
      ${optionalString bamOutput "| samtools view -b"} \
      > $out
  '';
  passthru.filetype = if bamOutput then filetype.bam { inherit ref; sorting = sort.none { }; } else filetype.sam { inherit ref; sorting = sort.name { }; };
  passthru.multicore = true;
}
