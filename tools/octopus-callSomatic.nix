{ bionix
, faidxAttrs ? {}
, indexAttrs ? {}
, flags ? ""}:

with bionix;
with lib;
with types;

{normal, tumours}:

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

  inputs = [normal] ++ tumours;
  getref = f: matchFiletype "octopus-callSomatic" { bam = {ref, ...}: ref; cram = {ref, ...}: ref;} f;
  refs = map getref inputs;
  ref = head refs;

in

assert (length (unique refs) == 1);


stage {
  name = "octopus-callSomatic";
  buildInputs = with pkgs; [ octopus-caller samtools ];
  buildCommand = ''
    ln -s ${ref} ref.fa
    ln -s ${samtools.faidx faidxAttrs ref} ref.fai
    ${concatMapStringsSep "\n" (i: ''
      ln -s ${i} $(basename ${i}).bam
      ln -s ${samtools.index indexAttrs i} $(basename ${i}).bai
    '') inputs}
    normal=$(samtools view -H ${normal} | awk -f ${smScript})
    octopus -R ref.fa -I *.bam -o $out \
      --threads=$NIX_BUILD_CORES \
      -N $normal \
      ${flags}
  '';
  passthru.filetype = filetype.vcf {ref = ref;};
  passthru.multicore = true;
}
