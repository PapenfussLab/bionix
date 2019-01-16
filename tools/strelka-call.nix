{ bionix
, indexAttrs ? {}
, bamIndexAttrs ? {}
, flags ? null
}:

inputs:

with bionix;
with lib;
with types;

let
  filename = path: last (splitString "/" path);
  getref = f: matchFiletype "strelka-call" { bam = x: x.ref; } f;
  refs = map getref inputs;
  ref = head refs;

  drv = bionix.strelka.call {inherit indexAttrs bamIndexAttrs flags;} inputs;

in

assert (length (unique refs) == 1);

stage {
  name = "strelka-call";
  buildInputs = with pkgs; [ strelka ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${bionix.samtools.faidx indexAttrs ref} ref.fa.fai
    ${concatMapStringsSep "\n" (p: "ln -s ${p} ${filename p}.bam") inputs}
    ${concatMapStringsSep "\n" (p: "ln -s ${bionix.samtools.index bamIndexAttrs p} ${filename p}.bai") inputs}

    configureStrelkaGermlineWorkflow.py \
    ${concatMapStringsSep " " (i: "--bam ${filename i}.bam") inputs} \
      --ref ref.fa \
      --runDir $TMPDIR

    ./runWorkflow.py \
      -m local \
      -j $NIX_BUILD_CORES 2>&1

    cp -r results $out
  '';
  passthru.variants = stage {
    name = "strelka-call-variants";
    buildCommand = ''
      ln -s ${drv}/variants/variants.vcf.gz $out
    '';
    passthru.filetype = filetype.gz (filetype.vcf {ref=ref;});
  };
}
