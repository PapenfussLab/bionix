{bionix, nixpkgs}:

with bionix;

{
  callVariants = callBionix ./gridss-callVariants.nix;
}
