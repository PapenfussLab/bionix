{bionix
,flags ? null}:

with bionix;
with lib;

{ref
,expr
,pos}:

stage {
  name = "inferCNV";
  buildInputs = [ bionix.infercnv.app ];
  buildCommand = ''
    inferCNV.R --output_dir $TMPDIR \
      ${optionalString (flags != null) flags} \
      --ref ${ref} \
      ${expr} \
      ${pos}
    mkdir $out
    cp -r $TMPDIR/* $out
  '';
}
