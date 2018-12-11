{ bionix
, nixpkgs
, nameSort ? false
, flags ? null
, outfmt ? null
}:

input:

with nixpkgs;
with lib;

let
  inherit (bionix.types) matchFiletype coordSort matchFileSorting;
in

assert (matchFiletype "samtools-sort" { bam = _: true; sam = _: true; cram = _: true; } input);

let
  outfmtR = if outfmt != null then outfmt input else input.filetype;
  outFmtFlags = matchFiletype "samtools-sort-outfmt" { bam = _: "-O BAM"; sam = _: "-O SAM"; cram = ref: "-O CRAM -T ${ref}"; } {filetype = outfmtR;};
  alreadySorted = matchFileSorting "samtools-sort" { name = _: nameSort; coord = _: !nameSort; none = _: false;} input;
in stdenv.mkDerivation {
  name = "samtools-sort";
  buildInputs = [ samtools ];
  buildCommand =
    if alreadySorted then
      "ln -s ${input} $out"
    else
      ''
        samtools sort -@ $NIX_BUILD_CORES \
          ${optionalString nameSort "-n"} \
          ${outFmtFlags} \
          ${optionalString (flags != null) flags} \
          ${input} > $out
      '';
  passthru.filetype = if nameSort then bionix.types.nameSort outfmtR else coordSort outfmtR;
}
