{bionix}:

with bionix;

{
  index = callBionixE ./kallisto-index.nix;

  /* Quantify reads against a transcriptome
  Type: { bias :: bool, bootstrapSamples :: int, seed :: int, plaintext :: bool, fusion :: bool, single :: bool, frStranded :: bool, rfStranded :: bool, fragmentLength :: Int, fragmentSD :: real } -> [fastq] -> kallisto
  */
  quant = callBionixE ./kallisto-quant.nix;
}
