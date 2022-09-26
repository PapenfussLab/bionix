{
  bionix,
  nameSort ? false,
  flags ? null,
  outfmt ? null,
  targets ? null,
}: input:
with bionix;
with lib;
with types;
assert (matchFiletype "samtools-view" {
    bam = _: true;
    sam = _: true;
    cram = _: true;
  }
  input); let
  handleTarget = x: let
    type = builtins.typeOf x;
    handler = handlers."${type}" or (builtins.throw "samtools-view:unhandled target type:${type}");
    handlers = {
      string = handleTarget [x];
      list = let
        file = pkgs.writeText "target.bed" (concatStringsSep "\n" x);
      in "--target-file ${file}";
      path = "--target-file ${x}";
      set = "--target-file ${x}";
    };
  in
    handler;

  outfmtR =
    if outfmt != null
    then
      (
        if builtins.typeOf outfmt == "string"
        then
          {
            "bam" = toBam;
            "cram" = toCram;
            "sam" = toSam;
          }
          ."${outfmt}"
        else outfmt
      )
      input
    else input.filetype;
  fa = ref: matchFiletype "samtools-view-ref" {fa = _: ref;} ref;
  outfmtFlags = matchFiletype "samtools-view-outfmt" {
    bam = _: "-O BAM";
    sam = _: "-O SAM";
    cram = x: "-O CRAM -T ${fa x.ref}";
  } {filetype = outfmtR;};
in
  stage {
    name = "samtools-view";
    buildInputs = with pkgs; [samtools];
    buildCommand = ''
      samtools view \
        ${optionalString (targets != null) (handleTarget targets)} \
        ${outfmtFlags} ${optionalString (flags != null) flags} ${input} > $out
    '';
    passthru.filetype = outfmtR;
  }
