{ bionix }:

with bionix;

{
  /* Compute coverage of a reference given an alignemnt
  Type: { ... } -> bam -> mosdepth
  */
  depth = callBionixE ./mosdepth-depth.nix;

  /* Plot sample coverages. Names are optional.
  Type: { ... } -> { inputs :: [mosdepth], names :: [string] } -> html
  */
  plot = callBionixE ./mosdepth-plot.nix;
}
