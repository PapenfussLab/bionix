{ bionix, nixpkgs }:

with nixpkgs;
with bionix;

{
  crumble = callPackage ./crumble-app.nix {};
}
