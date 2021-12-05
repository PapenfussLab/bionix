{ bionix }:

with bionix;
with lib;

rec {
  jar =
    let
      version = "2.13.0";
    in
    pkgs.fetchurl {
      url = "https://github.com/PapenfussLab/gridss/releases/download/v${version}/gridss-${version}-gridss-jar-with-dependencies.jar";
      sha256 = "sha256-fjBBI9Tz8vHtAdctAhUp8mM9oh6WHPUUavEdjixUVfU=";
    };

  /* Generate configuration file for GRIDSS. Takes attribute sets to GRIDSS ini style format.
    Type: genConfig :: attrSet -> ini file
  */
  genConfig = callBionix ./gridss-configFile.nix { };

  /* Invoke the callVariants tool
    Type: callVariants :: {blacklist :: drv = null, config :: ini = null, heapSize :: String = "31g", ...} -> [bam] -> variants
  */
  callVariants = callBionixE ./gridss-callVariants.nix;

  /* Invoke computeSamTags tool
    Type: computeSamTags :: {config :: ini = null, heapSize :: String = "1G", ...} -> bam -> bam
  */
  computeSamTags = callBionixE ./gridss-computeSamTags.nix;

  /* Invoke softClipsToSplitReads tool
    Type: softClipsToSplitReads :: {alignerStreaming :: Bool = false, config :: ini = null, heapSize :: String = "2G", ...} -> bam -> bam
  */
  softClipsToSplitReads = callBionixE ./gridss-softClipsToSplitReads.nix;

  /* Invoke collectMetrics tool
    Type: collectMetrics :: {thresholdCoverage :: Int = 10000, config :: ini = null, heapSize :: String = "1G", ...} -> bam -> metrics
  */
  collectMetrics = callBionixE ./gridss-collectMetrics.nix;

  /* Invoke extractSVReads tool
    Type: extractSVReads :: {unmappedReads :: Bool = false, minClipLength :: Int = 5, config :: ini = null, ...} -> bam -> bam
  */
  extractSVReads = callBionixE ./gridss-extractSVReads.nix;

  /* Invoke assembly tool
    Type: assemble :: {config :: ini = null, heapSize :: String = "31g", ...} -> [bam] -> bam
  */
  assemble = callBionixE ./gridss-assemble.nix;
  shardedAssemble = n: a: input:
    let
      assemblies = genList
        (i:
          bionix.gridss.assemble
            (a // {
              jobNodes = n;
              jobIndex = i;
            })
            input)
        n;
    in
    if n <= 1 then
      bionix.gridss.assemble a input
    else
      bionix.gridss.assemble (a // { workdirs = map (a: a.work) assemblies; })
        input;

  /* Invoke identifyVariants tool
    Type: identifyVariants :: {config :: ini = null, heapSize :: String = "4g", ...} -> [bam] -> VCF
  */
  identifyVariants = exec
    (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).identify);

  /* Invoke annotateVariants tool
    Type: annotateVariants :: {config :: ini = null, heapSize :: String = "4g", ...} -> [bam] -> VCF
  */
  annotateVariants = exec
    (attrs: input: ((callBionix ./gridss-variants.nix attrs) input).annotate);

  # As annotateVariants except include assembly in output
  annotateAndAssemble = exec (attrs: input:
    ((callBionix ./gridss-variants.nix attrs) input).annotateAndAssemble);

  /* Preprocess BAM files to extract SV reads
    Type: preprocessBam :: bam -> bam
  */
  preprocessBam = with samtools;
    flip pipe [
      (gridss.extractSVReads { })
      (sort { nameSort = true; })
      (gridss.computeSamTags { })
      (sort { })
    ];

  /* Call SVs: entire pipeline including preprocessing. It is recommended to use this function rather than the individual above tools.
    Type: [bam] -> GRIDSS result
  */
  call = inputs: gridss.annotateVariants { } (map gridss.preprocessBam inputs);

  # As call but include assemblies in output
  callAndAssemble = inputs:
    gridss.annotateAndAssemble { } (map gridss.preprocessBam inputs);
}
