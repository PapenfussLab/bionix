{nixpkgs, bionix}:

with nixpkgs;
with bionix;

{
  app = callPackage ./infercnv-app.nix {};
  infercnv = callBionix ./infercnv-infer.nix {};
}
