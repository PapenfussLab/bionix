{ bionix }:

with bionix;

{
  call = callBionixE ./octopus-call.nix;
  callSomatic = callBionixE ./octopus-callSomatic.nix;
}
