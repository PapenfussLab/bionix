{ bionix, nixpkgs }:

with nixpkgs;

{
  depth = attrs: callPackage ./mosdepth-depth.nix attrs;
  plot = attrs: callPackage ./mosdepth-plot.nix attrs;
}
