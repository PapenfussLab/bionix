{bionix, nixpkgs}:

with nixpkgs;
with bionix;

rec {
  jar = fetchurl {
    url = "https://github.com/PapenfussLab/gridss/releases/download/v2.0.0/gridss-2.0.0-gridss-jar-with-dependencies.jar";
    sha256 = "01srl3qvv060whqg1y1fpxjc5cwga5wscs1bmf1v3z87dignra7k";
  };
  callVariants = callBionix ./gridss-callVariants.nix;
  computeSamTags = callBionix ./gridss-computeSamTags.nix;
  softClipsToSplitReads = callBionix ./gridss-softClipsToSplitReads.nix;
  collectMetrics = callBionix ./gridss-collectMetrics.nix;
  extractSVReads = callBionix ./gridss-extractSVReads.nix;
  assemble = callBionix ./gridss-assemble.nix;
  identifyVariants = callBionix ./gridss-identifyVariants.nix;
  annotateVariants = callBionix ./gridss-annotateVariants.nix;
  preprocessBam = input: with samtools; markdup {} (sort {} (fixmate {mateScore = true;} (softClipsToSplitReads {} (computeSamTags {} (sort {nameSort = true;} input)))));
  call = inputs: annotateVariants {} (map preprocessBam inputs);
}
