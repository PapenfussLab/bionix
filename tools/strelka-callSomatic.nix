{ bionix
, indexAttrs ? { }
, bamIndexAttrs ? { }
, flags ? null
}:

{ normal, tumour }:

with bionix;
with lib;
with types;

let
  filename = path: last (splitString "/" path);
  getref = matchFiletype "strelka-callSomatic" { bam = x: x.ref; };
  inputs = [ normal tumour ];
  refs = map getref inputs;
  ref = head refs;

in

assert (length (unique refs) == 1);

let

  out = stage {
    name = "strelka-callSomatic";
    buildInputs = with pkgs; [ strelka ];
    outputs = [ "out" "indels" "snvs" ];
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
        sed -i '/^##cmd/d' $g
      done
      mv somatic.indels.vcf $indels
      mv somatic.snvs.vcf $snvs

      ln -s $snvs $out
    '';
    passthru.multicore = true;
    passthru.filetype = types.filetype.vcf { inherit ref; };
  };

in
out
