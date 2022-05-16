{bionix}:
with bionix;
with lib.types;
  input: let
    bamft = matchFiletype "unquip" {qp = x: x;} input;
    ref = matchFiletype' "unquip-sub" {bam = f: f.ref;} bamft;
  in
    stage {
      name = "unquip";
      buildInputs = [quip.app];
      buildCommand = ''
        ln -s ${input} input.bam.qp
        unquip -r ${ref} input.bam.qp
        cp input.bam $out
      '';
      stripStorePaths = false;
      passthru.filetype = bamft;
    }
