{ bionix
, fast ? false
, very-fast ? false
, max-genotypes ? null
, targets ? null
, faidxAttrs ? { }
, indexAttrs ? { }
, flags ? ""
}:

assert !fast || !very-fast;
assert max-genotypes == null || max-genotypes > 0;

with bionix;
with lib;
with types;

inputs:

let
  getref = matchFiletype "octopus-call" { bam = { ref, ... }: ref; cram = { ref, ... }: ref; };
  refs = map getref inputs;
  ref = head refs;


  handleTarget = x:
    let
      type = builtins.typeOf x;
      handler = handlers."${type}" or (builtins.throw "octopus-callSomatic:unhandled target type:${type}");
      handlers = {
        string = "-T '${x}'";
        list =
          let file = pkgs.writeText "regions.txt" (concatStringsSep "\n" x);
          in "-t ${file}";
        path = "-t ${x}";
        set = "-t ${x}";
      };
    in
    handler;


in

assert (length (unique refs) == 1);

stage {
  name = "octopus-call";
  buildInputs = with pkgs; [ octopus-caller ];
  outputs = [ "out" "evidence" ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fai
    ${concatMapStringsSep "\n" (i: ''
      ln -s ${i} $(basename ${i}).bam
      ln -s ${samtools.index indexAttrs i} $(basename ${i}).bai
    '') inputs}
    ${optionalString (length inputs > 1) "mkdir $evidence"}
    octopus -R ref.fa -I *.bam -o $out \
      --bamout $evidence \
      --threads=$NIX_BUILD_CORES \
      ${optionalString fast "--fast"} \
      ${optionalString very-fast "--very-fast"} \
      ${optionalString (max-genotypes != null) "--max-genotypes ${toString max-genotypes}"} \
      ${optionalString (targets != null) (handleTarget targets)} \
      ${flags}

    # Strip out octopus ARGV
    sed -i '/^##octopus=/d' $out
  '';
  passthru.filetype = filetype.vcf { inherit ref; };
  passthru.multicore = true;
  stripStorePaths = false;
}
