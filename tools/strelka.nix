{ bionix }:

with bionix;

{
  callSomatic = callBionixE ./strelka-callSomatic.nix;
  call = callBionixE ./strelka-call.nix;
}
