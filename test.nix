(import ./modules.nix {
  configuration = { bwa.align.ref = ./examples/ref.fa; };
  pkgs = import <nixpkgs> { };
}).bionix.bwa.align {
  input1 = ./examples/sample1-1.fq;
  input2 = ./examples/sample1-2.fq;
}
