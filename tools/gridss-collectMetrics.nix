{ bionix
, nixpkgs
, thresholdCoverage ? 10000
, flags ? null
}:

with nixpkgs;
with lib;
with bionix.types;

input:

let
  ref = matchFiletype "gridss-collectMetrics" { bam = x: x.ref; } input;
in


stdenv.mkDerivation rec {
  name = "gridss-collectMetrics";
  buildInputs = [ jre R ];
  buildCommand = ''
    mkdir $out
    java -Xmx1G -cp ${bionix.gridss.jar} \
			gridss.analysis.CollectGridssMetrics \
			I=${input}\
			O=$out/input \
      AS=true \
			THRESHOLD_COVERAGE=${toString thresholdCoverage}
  '';
}
