{ bionix }:

with bionix;

{
  /* Call structural variants
    Type: { ... } -> { normals :: [bam], tumour :: bam } -> manta
  */
  call = callBionixE ./manta-call.nix;
}
