{bionix
,nixpkgs
,flags ? null}:

with nixpkgs;
with lib;

{ref
,expr
,pos}:

stdenv.mkDerivation {
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
