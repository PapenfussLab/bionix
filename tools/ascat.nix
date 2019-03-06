{bionix}:

with bionix;

{
  app = pkgs.callPackage ./ascat-app.nix {};

  /* Generate GC correction file for ascatNGS.
  Type: {ref :: fasta, chrPrefix :: string, ...} -> (snps :: VCF) -> (gc :: GC)
  */
  gccorrect = callBionixE ./ascat-gccorrect.nix;

  /* Call CNVs using ascatNGS. Gender is a string as per ascatNGS docs (e.g., "XX").
  Type: {ref :: fasta, gc :: GC, ...} -> {tumour :: bam, normal :: bam, gender :: string} -> CNVs
  */
  callCNV = callBionixE ./ascat-callCNV.nix;
}
