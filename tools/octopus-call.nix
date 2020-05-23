{ bionix
, fast ? false
, very-fast ? false
, max-genotypes ? null
, targets ? null
, faidxAttrs ? {}
, indexAttrs ? {}
, flags ? ""}:

assert !fast || !very-fast;
assert max-genotypes == null || max-genotypes > 0;

with bionix;
with lib;
with types;

inputs:

let
  getref = f: matchFiletype "octopus-call" { bam = {ref, ...}: ref; cram = {ref, ...}: ref;} f;
  refs = map getref inputs;
  ref = head refs;

in

assert (length (unique refs) == 1);

stage {
  name = "octopus-call";
  buildInputs = with pkgs; [ octopus-caller ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fai
    ${concatMapStringsSep "\n" (i: ''
      ln -s ${i} $(basename ${i}).bam
      ln -s ${samtools.index indexAttrs i} $(basename ${i}).bai
    '') inputs}
    octopus -R ref.fa -I *.bam -o $out \
      --threads=$NIX_BUILD_CORES \
      ${optionalString fast "--fast"} \
      ${optionalString very-fast "--very-fast"} \
      ${optionalString (max-genotypes != null) "--max-genotypes ${toString max-genotypes}"} \
      ${optionalString (targets != null) (if builtins.typeOf targets == "list" then "-T ${concatStringsSep "," targets}" else "-t ${targets}")} \
      ${flags}
  '';
  passthru.filetype = filetype.vcf {ref = ref;};
  passthru.multicore = true;
}
