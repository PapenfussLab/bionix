{bionix}:

with bionix;

rec {
  jar = pkgs.fetchurl {
    url = "https://github.com/PapenfussLab/gridss/releases/download/v2.0.0/gridss-2.0.0-gridss-jar-with-dependencies.jar";
    sha256 = "01srl3qvv060whqg1y1fpxjc5cwga5wscs1bmf1v3z87dignra7k";
  };
  genConfig = callBionixE ./gridss-configFile.nix {};
  callVariants = callBionixE ./gridss-callVariants.nix;
  computeSamTags = callBionixE ./gridss-computeSamTags.nix;
  softClipsToSplitReads = callBionixE ./gridss-softClipsToSplitReads.nix;
  collectMetrics = callBionixE ./gridss-collectMetrics.nix;
  extractSVReads = callBionixE ./gridss-extractSVReads.nix;
  assemble = callBionixE ./gridss-assemble.nix;
  identifyVariants = exec (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).identify);
  annotateVariants = exec (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).annotate);
  preprocessBam = with samtools;
    pipe [
      (gridss.extractSVReads {})
      (sort {nameSort = true;})
      (gridss.computeSamTags {})
      (gridss.softClipsToSplitReads {})
      (sort {})
    ];
  call = inputs: gridss.annotateVariants {} (map gridss.preprocessBam inputs);
}
