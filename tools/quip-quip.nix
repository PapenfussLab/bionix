{bionix}:
with bionix;
with lib.types;
  input: let
    ref = matchFiletype "quip" {bam = f: f.ref;} input;
  in
    stage {
      name = "quip";
      buildInputs = [quip.app];
      buildCommand = ''
        ln -s ${input} input.bam
        quip -r ${ref} input.bam
        cp input.bam.qp $out
      '';
      stripStorePaths = false;
      passthru.filetype = filetype.qp input.filetype;
    }
