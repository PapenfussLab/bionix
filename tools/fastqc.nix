{ bionix }:

with bionix;

{
  fastqc = pkgs.callPackage ./fastqc-app.nix {};

  /* QC check
  Type: check :: {...} -> input :: fastq -> report
  */
  check = callBionixE ./fastqc-check.nix;
}
