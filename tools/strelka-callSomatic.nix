{ bionix
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

{normal, tumour}:

with bionix;
with lib;
with types;

let
  filename = path: last (splitString "/" path);
  getref = f: matchFiletype "strelka-callSomatic" { bam = x: x.ref; } f;
  inputs = [ normal tumour ];
  refs = map getref inputs;
  ref = head refs;

in

assert (length (unique refs) == 1);

stage {
  name = "strelka-callSomatic";
  buildInputs = with pkgs; [ strelka ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx indexAttrs ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bionix.samtools.index bamIndexAttrs p} ${filename p}.bai") inputs}

    configureStrelkaSomaticWorkflow.py \
      --normalBam ${filename normal}.bam \
      --tumourBam ${filename tumour}.bam \
      --ref ref.fa \
      --runDir $TMPDIR

    ./runWorkflow.py \
      -m local \
      -j $NIX_BUILD_CORES

    # Strelka writes runtime stats and timestamps;
    # both have to be stripped to provide determinism
    cd results/variants
    rm *.tbi
    for f in *.vcf.gz; do
      gunzip $f
      g=$(basename $f .gz)
      sed -i '/^##fileDate/d' $g
      sed -i '/^##startTime/d' $g
    done
    mkdir $out
    cp -r * $out
  '';
  passthru.multicore = true;
}
