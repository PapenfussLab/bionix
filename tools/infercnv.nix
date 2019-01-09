{bionix}:

with bionix;

{
  app = callPackage ./infercnv-app.nix {};
  infercnv = callBionixE ./infercnv-infer.nix {};
}
