{ bionix }:

with bionix;
with types;

{
  /* Calls somatic variants
    Type: callSomatic :: {...} -> {tumour, normal} -> somatic results
  */
  callSomatic = callBionixE ./strelka-callSomatic.nix;
  /* Calls variants
    Type: call :: {...} -> [input] -> results
  */
  call = callBionixE ./strelka-call.nix;
}
