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

{ normal, tumours }:

let
  smScript = pkgs.writeText "smScript.awk" ''
    /^@RG/{
      for(i = 1; i <= NF; i++) {
        n=split($i, fields, ":")
        if(n == 2 && fields[1] == "SM"){
          print fields[2]
          exit
        }
      }
    }
  '';

  inputs = [ normal ] ++ tumours;
  getref = f: matchFiletype "octopus-callSomatic" { bam = { ref, ... }: ref; cram = { ref, ... }: ref; } f;
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
      };
    in
    handler;

in

assert (length (unique refs) == 1);


stage {
  name = "octopus-callSomatic";
  buildInputs = with pkgs; [ octopus-caller samtools ];
  outputs = [ "out" "evidence" ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fai
    ${concatMapStringsSep "\n" (i: ''
      ln -s ${i} $(basename ${i}).bam
      ln -s ${samtools.index indexAttrs i} $(basename ${i}).bai
    '') inputs}
    normal=$(samtools view -H ${normal} | awk -f ${smScript})
    mkdir $evidence
    octopus -R ref.fa -I *.bam -o $out \
      --bamout $evidence \
      --threads=$NIX_BUILD_CORES \
      ${optionalString fast "--fast"} \
      ${optionalString very-fast "--very-fast"} \
      ${optionalString (max-genotypes != null) "--max-genotypes ${toString max-genotypes}"} \
      ${optionalString (targets != null) (handleTarget targets)} \
      -N $normal \
      ${flags}
  '';
  passthru.filetype = filetype.vcf { ref = ref; };
  passthru.multicore = true;
}
