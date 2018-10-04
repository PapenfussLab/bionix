{ bionix
, nixpkgs
, ref
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

{normal, tumour}:

with nixpkgs;
with lib;

let
  filename = path: last (splitString "/" path);
  inputs = [ normal tumour ];

in stdenv.mkDerivation {
  name = "strelka";
  buildInputs = [ strelka ];
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
