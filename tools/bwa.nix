{ bionix, nixpkgs }:

with nixpkgs;

{
  align = attrs: callPackage ./bwa-mem.nix attrs;
  index = attrs: callPackage ./bwa-index.nix attrs;
}
