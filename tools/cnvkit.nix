{bionix}:

with bionix;

{
  /* Call CNVs
  Type: callCNV :: {targets :: target file, annotations :: annotation file, ...} -> {normals :: [bam], tumours :: [bam]} -> CNVs
  */
  callCNV = callBionixE ./cnvkit-batch.nix;

  /* Scatter plot from CNV calls
  Type: scatterPlot :: {} -> CNVs -> PDF
  */
  scatterPlot = callBionixE ./cnvkit-scatter.nix;
}
