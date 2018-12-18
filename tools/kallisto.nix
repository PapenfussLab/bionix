{bionix, nixpkgs}:

with bionix;

{
  index = callBionixE ./kallisto-index.nix;
  quant = callBionixE ./kallisto-quant.nix;
}
