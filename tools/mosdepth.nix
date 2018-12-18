{ bionix, nixpkgs }:

with bionix;

{
  depth = callBionixE ./mosdepth-depth.nix;
  plot = callBionixE ./mosdepth-plot.nix;
}
