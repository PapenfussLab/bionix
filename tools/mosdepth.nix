{ bionix, nixpkgs }:

with bionix;

{
  depth = callBionix ./mosdepth-depth.nix;
  plot = callBionix ./mosdepth-plot.nix;
}
