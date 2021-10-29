{ bionix }:

with bionix;
with lib.types;

{
  regex = n: input: outputDrvs (callBionixE ./shard-regex.nix { inherit n; } input);
  fastQPair = n: { input1, input2 }: lib.zipListsWith (i: j: { input1 = i; input2 = j; }) (lib.shard.regex n input1) (lib.shard.regex n input2);
}
