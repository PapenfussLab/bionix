{ bionix
, nixpkgs
, thresholdCoverage ? 10000
, flags ? null
, config ? null
}:

with nixpkgs;
with lib;
with bionix.types;

input:

let
  ref = matchFiletype "gridss-collectGridssMetrics" { bam = x: x.ref; } input;
  sorted = matchFileSorting "gridss-collectGridssMetrics" { name = _: true; } input;
in


stdenv.mkDerivation rec {
  name = "gridss-collectGridssMetrics";
  buildInputs = [ jre ];
  buildCommand = ''
    mkdir $out
    ln -s ${input} input.bam
    java -Xmx1G -cp ${bionix.gridss.jar} \
			gridss.analysis.CollectGridssMetrics \
			${optionalString sorted "ASSUME_SORTED=true"} \
      ${optionalString (config != null) ("CONFIGURATION_FILE=" + bionix.gridss.gridssConfig config)} \
			I=input.bam \
			O=$out \
			THRESHOLD_COVERAGE=${toString thresholdCoverage}
    '';
}
