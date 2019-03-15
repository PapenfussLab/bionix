{bionix
, indexFlags ? {}
, bias ? false
, bootstrapSamples ? 0
, seed ? 42
, plaintext ? false
, fusion ? false
, single ? false
, frStranded ? false
, rfStranded ? false
, fragmentLength ? null
, fragmentSD ? null
, ref}:

with bionix;
with lib;

assert (!single || (fragmentLength != null && fragmentSD != null));

inputs:

let
  inherit (bionix.types) matchFiletype';
  isFastQ = matchFiletype' "kallisto-quant" {fq = _: true; gz = isFastQ; };
in

assert (all (x: isFastQ (x.filetype)) inputs);

stage {
  name = "kallisto-quant";
  buildInputs = with pkgs; [ kallisto ];
  buildCommand = ''
    mkdir $out
    kallisto quant \
      -i ${bionix.kallisto.index indexFlags ref} \
      -o $out \
      ${optionalString bias "--bias"} \
      ${optionalString (bootstrapSamples > 0) "-b ${toString bootstrapSamples} --seed=${toString seed}"} \
      ${optionalString plaintext "--plaintext"} \
      ${optionalString fusion "--fusion"} \
      ${optionalString single "--single -l ${toString fragmentLength} -s ${toString fragmentSD}"} \
      ${optionalString frStranded "--fr-stranded"} \
      ${optionalString rfStranded "--rf-stranded"} \
      -t $NIX_BUILD_CORES \
      ${concatStringsSep " " inputs}
  '';
  passthru.multicore = true;
}
