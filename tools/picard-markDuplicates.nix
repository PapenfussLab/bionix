{ bionix
, flags ? null
} :

inputBam :

with bionix;
with lib;
with types;

assert (matchFiletype "picard-markDuplicates" { bam = _: true; } inputBam);

# Not sure what to do with sorting: behavior varies based on sortedness of input!
# See: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
# assert (matchFileSorting "picard-markDuplicates" { coord = _: true; } inputBam);

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
