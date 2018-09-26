{bionix, nixpkgs}:

with nixpkgs;

{
  callVariants = attrs: callPackage ./gridss-callVariants.nix attrs;
}
