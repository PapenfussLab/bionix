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
  name = "strelka";
  buildInputs = with pkgs; [ strelka gzip ];
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

    cp -r results $out
  '';
}
