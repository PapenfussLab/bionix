{ bionix, nixpkgs }:

with nixpkgs;

{
  call = attrs: callPackage ./strelka-call.nix attrs;
}
