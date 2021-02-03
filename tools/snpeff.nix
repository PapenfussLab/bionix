{bionix}:

with bionix;

{
  /* Annotate variants with SNPEff database. Some annotation DBs available in bionix.ref (e.g., bionix.ref.grch38.snpeff.db).
  Type: { db, ... } -> vcf -> vcf
  */
  annotate = callBionixE ./snpeff-annotate.nix;
  
  /* Annotate variants with dbNSFP database. Some dbNSFP annotation DBs available in bionix.ref (e.g., bionix.ref.grch38.snpeff.dbNSFP).
  Type: { dbnsfp, ... } -> vcf -> vcf
  */
  dbnsfp = callBionixE ./snpeff-dbnsfp.nix;
}
