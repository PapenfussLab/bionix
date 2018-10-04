{bionix, nixpkgs}:

with bionix;

{
  callVariants = callBiolnix ./gridss-callVariants.nix;
}
