{ bionix }:

with bionix;

{
  fastqc = pkgs.callPackage ./fastqc-app.nix {};
  check = callBionixE ./fastqc-check.nix;
}
