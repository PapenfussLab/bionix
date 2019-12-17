{ bionix, indexAttrs ? {}, faidxAttrs ? {}, flags ? "" }:

{ normals ? [], tumour ? null }:

with bionix;
with lib;
with types;

let
  getref = matchFiletype "manta-call" { bam = x: x.ref; };
  refs = map getref normals ++ optionals (tumour != null) [(getref tumour)];
  ref = head refs;

  renameAndIndex = f: stage {
    name = "rename";
    buildCommand = ''
      mkdir $out
      ln -s ${f} $out/sample.bam
      ln -s ${samtools.index indexAttrs f} $out/sample.bam.bai
    '';
  };
in

assert (length (unique refs) == 1);

stage {
  name = "manta-call";
  buildInputs = with pkgs; [ manta strace ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx faidxAttrs ref} ref.fa.fai
    configManta.py ${optionalString (normals != null) (concatMapStringsSep " " (n: "--bam=${renameAndIndex n}/sample.bam") normals)} \
      ${optionalString (tumour != null) "--tumourBam=${renameAndIndex tumour}/sample.bam"} \
      --runDir=$TMPDIR \
      --referenceFasta=ref.fa \
      ${flags}
    ./runWorkflow.py -j $NIX_BUILD_CORES
    cp -r results $out
  '';
  passthru.multicore = true;
}
