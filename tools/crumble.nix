{ bionix, nixpkgs }:

with nixpkgs;
with bionix;

{
  crumble = callPackage ./crumble-app.nix {};
  toCram = callBionix ./crumble-toCram.nix;
}
