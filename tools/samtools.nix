{ bionix, nixpkgs }:

with nixpkgs;

{
  index = attrs: callPackage ./samtools-index.nix attrs;
  sort = attrs: callPackage ./samtools-sort.nix attrs;
  faidx = attrs: callPackage ./samtools-faidx.nix attrs;
}
