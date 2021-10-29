{ bionix
, flags ? null
, ploidy
}:

inputs:

with bionix;
with lib;
with types;

let
  config = pkgs.writeText "pool.txt" (concatMapStringsSep "\n" (x: "${x}\t${toString ploidy}\t1") (attrValues inputs));
  getref = matchFiletype "SNVer-call" { bam = { ref, ... }: ref; };
  refs = mapAttrsToList (_: getref) inputs;
  ref = head refs;
in

assert (length (unique refs) == 1);

stage {
  name = "SNVerPool";
  buildInputs = [ snver.app ];
  outputs = [ "out" "log" "raw" "filter" "indelfilter" "indelraw" ];
  buildCommand = ''
    SNVerPool -i / -c ${config} -r ${ref} -o snver

    for f in *.vcf ; do
      ${concatStringsSep "\n" (mapAttrsToList (x: y: "sed -i 's|//${y}|${x}|' $f") inputs)}
    done

    mkdir $out
    cp snver.failed.log $log
    ln -s $log $out/snver.failed.log
    cp snver.filter.vcf $filter
    ln -s $filter $out/snver.filter.vcf
    cp snver.indel.filter.vcf $indelfilter
    ln -s $indelfilter $out/snver.indel.filter.vcf
    cp snver.indel.raw.vcf $indelraw
    ln -s $indelraw $out/snver.indel.raw.vcf
    cp snver.raw.vcf $raw
    ln -s $raw $out/snver.raw.vcf
  '';
}
