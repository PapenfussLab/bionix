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
  inherit (bionix.types) matchFiletype coordSort;
in

assert (matchFiletype "samtools-sort" { bam = _: true; sam = _: true; cram = _: true; } input);

let
  outfmtR = if outfmt != null then outfmt input else input.filetype;
  outFmtFlags = matchFiletype "samtools-sort-outfmt" { bam = _: "-O BAM"; sam = _: "-O SAM"; cram = ref: "-O CRAM -T ${ref}"; } {filetype = outfmtR;};
in stdenv.mkDerivation {
  name = "samtools-sort";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools sort -@ $NIX_BUILD_CORES ${optionalString nameSort "-n"} ${outFmtFlags} ${optionalString (flags != null) flags} ${input} > $out
  '';
  passthru.filetype = if nameSort then bionix.types.nameSort outfmtR else coordSort outfmtR;
}
