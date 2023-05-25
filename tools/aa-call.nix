{
  bionix,
  seeds,
  license,
  flags ? "",
  downsample ? 0.1,
  ...
}: input:
with bionix; let
  indexedBam = linkOutputs {
    "input.bam" = input;
    "input.bam.bai" = samtools.index {} input;
  };
in
  stage {
    name = "aa-call";
    MOSEKLM_LICENSE_FILE = license;
    buildInputs = [bionix.ampliconarchitect.app];
    buildCommand = ''
      mkdir $out
      export AA_DATA_REPO=$TMPDIR
      tar -xzf ${self.aa.ref}
      AmpliconArchitect.py \
        --bam ${indexedBam}/input.bam \
        --bed ${seeds} \
        --ref GRCh38 \
        --out $out/out \
        --downsample ${toString downsample} \
        $flags
    '';
  }
