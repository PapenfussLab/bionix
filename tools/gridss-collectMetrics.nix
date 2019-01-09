{ bionix
, thresholdCoverage ? 10000
, flags ? null
, config ? null
, heapSize ? "1G"
}:

with bionix;
with lib;
with types;

input:

let
  ref = matchFiletype "gridss-collectMetrics" { bam = x: x.ref; } input;
in


stage rec {
  name = "gridss-collectMetrics";
  buildInputs = with pkgs; [ jre R ];
  buildCommand = ''
    mkdir $out
    java -Xmx${heapSize} -cp ${bionix.gridss.jar} \
			gridss.analysis.CollectGridssMetrics \
      ${optionalString (config != null) ("OPTIONS_FILE=" + bionix.gridss.gridssConfig config)} \
			I=${input}\
			O=$out/input \
      AS=true \
			THRESHOLD_COVERAGE=${toString thresholdCoverage}
  '';
}
