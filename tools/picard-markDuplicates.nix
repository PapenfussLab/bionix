{ bionix
, flags ? null
}:

inputBam:

with bionix;
with lib;
with types;

assert (matchFiletype "picard-markDuplicates" { bam = _: true; } inputBam);
assert (matchFileSorting
  "picard-markDuplicates"
  { coord = _: true; name = _: true; }
  inputBam);

# Note that picard markDuplicates has different behaviour depending on whether the input
# is name-sorted or coordinate-sorted.

stage {
  name = "picard-markDuplicates";
  buildInputs = with pkgs;
    [ picard-tools ];
  outputs = [ "out" "metrics" ];
  buildCommand = ''
    picard MarkDuplicates \
        I=${inputBam} \
        O=$out \
        M=$metrics \
        ${optionalString (flags != null) flags}
  '';
  passthru.filetype = inputBam.filetype;
}
