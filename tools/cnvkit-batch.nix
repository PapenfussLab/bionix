{ bionix
, targets ? null
, annotations ? null
, flags ? null
, indexAttrs ? { }
}:

{ normals ? [ ], tumours }:

with bionix;
with lib;
with types;

let
  getref = matchFiletype "cnvkit-batch" { bam = { ref, ... }: ref; };
  refs = map getref normals ++ map getref tumours;
  ref = head refs;
  sorted = matchFileSorting "cnvkit-batch" { coord = _: true; };
in

assert (length (unique refs) == 1);
assert (all sorted (normals ++ tumours));

stage {
  name = "cnvkit";
  buildInputs = with pkgs; [ python3Packages.cnvkit ];
  outputs = [ "out" ] ++ builtins.genList (x: "out${toString (x + 1)}") (length tumours);
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

    # Copy individual tumour files
    mkdir $out
    cnt=1
    for f in ${concatStringsSep " " tumours} ; do
      output="out$cnt"
      mkdir ''${!output}
      for g in $(basename $f)*.{cnr,cnn,cns} ; do
        cp $g ''${!output}/sample-''${g#*-}
      done
      cnt=$((cnt+1))
      ln -s ''${!output} $out/$output
    done
  '';
  passthru.multicore = true;
}
