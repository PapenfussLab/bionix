{ bionix
, nixpkgs
, nameSort ? false
, flags ? null
, outfmt ? null
}:

input:

with nixpkgs;
with lib;
with bionix.types;

assert (matchFiletype "samtools-view" { bam = _: true; sam = _: true; cram = _: true; } input);

let
  outfmtR = if outfmt != null then outfmt input else input.filetype;
  fa = ref: matchFiletype "samtools-view-ref" { fa = _: ref; } ref;
  outfmtFlags = matchFiletype "samtools-view-outfmt" { bam = _: "-O BAM"; sam = _: "-O SAM"; cram = x: "-O CRAM -T ${fa x.ref}"; } {filetype = outfmtR;};
in stdenv.mkDerivation {
  name = "samtools-view";
  buildInputs = [ samtools ];
  buildCommand = ''
    samtools view ${outfmtFlags} ${optionalString (flags != null) flags} ${input} > $out
  '';
  passthru.filetype = outfmtR;
}
