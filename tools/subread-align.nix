{ bionix
, ref
, flags ? null
, indexAttrs ? { }
, RG ? { }
, type ? 1
}:

{ input1
, input2 ? null
}:

with bionix;
with lib;
with types;
with compression;

let
  fa = f: matchFiletype "subread-ref" { fa = _: f; } f;
  fq = f: matchFiletype "subread-input" { fq = _: f; gz = matchFiletype' "subread-input" { fq = _: f; }; } f;

in
stage {
  name = "subread-align";
  buildInputs = with pkgs; [ subread ];
  buildCommand = ''
    subread-align -i ${bionix.subread.index indexAttrs ref}/ref ${optionalString (flags != null) flags} -T $NIX_BUILD_CORES \
      -t ${toString type} \
      -r ${fq input1} \
      ${optionalString (input2 != null) "-R ${fq input2}"} \
      ${optionalString (RG ? ID) ''
        --rg-id ${RG.ID} ${concatMapAttrsStringsSep " " (k: v: "--rg ${k}:${v}") (filterAttrs (k: _: k != "ID") RG)} \
      ''} \
      -o $out
  '';
  passthru.filetype = filetype.bam { inherit ref; sorting = sort.none { }; };
  passthru.multicore = true;
}
