{bionix}:

with bionix;

{
  /* Split aligned reads based on primary and secondary alignments
  Type: { ... } -> { primary :: bam, secondary :: bam } -> { primary_specific :: bam, primary_multi :: bam, secondary_specific :: bam, secondary_multi :: bam, unassigned :: bam, unresolved :: bam }
  */
  allocate = callBionixE ./xenomapper-allocate.nix;
}