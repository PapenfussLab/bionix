{bionix, nixpkgs}:

with nixpkgs;
with bionix;

rec {
  jar = fetchurl {
    url = "https://github.com/PapenfussLab/gridss/releases/download/v2.0.0/gridss-2.0.0-gridss-jar-with-dependencies.jar";
    sha256 = "01srl3qvv060whqg1y1fpxjc5cwga5wscs1bmf1v3z87dignra7k";
  };
  gridssConfig = callBionixE ./gridss-configFile.nix {};
  callVariants = callBionixE ./gridss-callVariants.nix;
  computeSamTags = callBionixE ./gridss-computeSamTags.nix;
  softClipsToSplitReads = callBionixE ./gridss-softClipsToSplitReads.nix;
  collectMetrics = callBionixE ./gridss-collectMetrics.nix;
  extractSVReads = callBionixE ./gridss-extractSVReads.nix;
  assemble = callBionixE ./gridss-assemble.nix;
  identifyVariants = exec (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).identify);
  annotateVariants = exec (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).annotate);
  preprocessBam = input: with samtools; sort {} (gridss.softClipsToSplitReads {} (gridss.computeSamTags {} (sort {nameSort = true;} (gridss.extractSVReads {} (markdup {} (sort {} (fixmate {mateScore = true;} (sort {nameSort = true;} input))))))));
  call = inputs: bionix.gridss.annotateVariants {} (map gridss.preprocessBam inputs);
}
