{bionix
,targets ? null
,annotations ? null
,flags ? null
,indexAttrs ? {}}:

{normals ? [], tumours}:

with bionix;
with lib;
with types;

let
  getref = f: matchFiletype "cnvkit-batch" { bam = {ref, ...}: ref; } f;
  refs = map getref normals ++ map getref tumours;
  ref = head refs;
  sorted = matchFileSorting "cnvkit-batch" { coord = _: true; };
in

assert (length (unique refs) == 1);
assert (all sorted (normals ++ tumours));

stage {
  name = "cnvkit";
  buildInputs = [ cnvkit.app ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx indexAttrs ref} ref.fa.fai
    cnvkit.py batch ${concatStringsSep " " tumours} \
      ${optionalString (normals != []) ("-n " + concatStringsSep " " normals)} \
      ${optionalString (annotations != null) "--annotate ${annotations}"} \
      ${if targets != null then "--targets ${targets}" else "-m wgs"} \
      -f ref.fa \
      -p $NIX_BUILD_CORES \
      -d $TMPDIR \
      ${optionalString (flags != null) flags}
    mkdir $out
    cp * $out
  '';
}
