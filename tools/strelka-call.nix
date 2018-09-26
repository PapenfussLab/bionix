{ stdenv
, callPackage
, lib
, strelka
, ref
, index ? callPackage ./samtools-faidx.nix {}
, bamIndex ? callPackage ./samtools-index.nix {}
, flags ? null
}:

{normal, tumour}:

with lib;

let
  filename = path: last (splitString "/" path);
  inputs = [ normal tumour ];

in stdenv.mkDerivation {
  name = "strelka";
  buildInputs = [ strelka ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${index ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bamIndex p} ${filename p}.bai") inputs}

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
