{pkgs ? import <nixpkgs> {}
,bionix ? import <bionix> {}}:

with pkgs;
with lib;
with bionix;

let
  inherit (types) filetype tagFiletype;

  fetchLocal = path: stdenv.mkDerivation {
    name = baseNameOf path;
    buildCommand = "ln -s ${path} $out";
  };
  tagfq = path: tagFiletype (filetype.fq {}) (fetchLocal path);
  tagfa = path: tagFiletype (filetype.fa {}) (fetchLocal path);

  ref = tagfa ./example/ref.fa;
  alignWithRG = rg: bwa.align { inherit ref; flags = "-R'@RG\\tID:${rg}\\tSM:${rg}'";};

  samples = [ {name = "mysample1"; files = {input1 = tagfq ./example/sample1-1.fq; input2 = tagfq ./example/sample1-2.fq;};}
              {name = "mysample2"; files = {input1 = tagfq ./example/sample2-1.fq; input2 = tagfq ./example/sample2-1.fq;};} ];

  alignments = map (i: samtools.sort {} (alignWithRG i.name i.files)) samples;
  variants = platypus.call {} alignments;

in stdenv.mkDerivation {
  name = "myproject";
  buildCommand = ''
    mkdir $out
    ln -s ${variants} $out/platypus.vcf
    mkdir $out/alignments
    ${concatStringsSep "\n" (zipListsWith (s: a: "ln -s ${a} $out/alignments/${s.name}.bam") samples alignments)}
  '';
}
