{bionix, nixpkgs}:

with bionix;

{
  index = callBionix ./kallisto-index.nix;
  quant = callBionix ./kallisto-quant.nix;
}
