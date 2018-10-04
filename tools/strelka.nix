{ bionix, nixpkgs }:

with bionix;

{
  call = callBionix ./strelka-call.nix;
}
