{ bionix, nixpkgs }:

with bionix;

{
  call = callBionixE ./strelka-call.nix;
}
