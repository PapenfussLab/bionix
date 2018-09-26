{ bionix, nixpkgs }:

with nixpkgs;

{
  faidx = attrs: callPackage ./samtools-faidx.nix attrs;
  flagstat = attrs: callPackage ./samtools-flagstat.nix attrs;
  index = attrs: callPackage ./samtools-index.nix attrs;
  sort = attrs: callPackage ./samtools-sort.nix attrs;
}
