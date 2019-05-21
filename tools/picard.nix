{ bionix }:

with bionix;

{
    markDuplicates = callBionixE ./picard-markDuplicates.nix;
}