{ bionix
, flags ? null
} :

inputBam :

with bionix;
with lib;
with types;

assert (matchFiletype "picard-markDuplicates" { bam = _: true; } input);
assert (matchFileSorting "picard-markDuplicates" { coord = _: true; } input);

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
}