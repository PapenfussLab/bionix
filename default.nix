{nixpkgs ? import <nixpkgs> {}}:

let
  bionix = nixpkgs.lib.makeExtensible (self:
    let callBionix = file: import file { bionix = self; nixpkgs = nixpkgs; };
    in with self; {
      bwa = callBionix ./tools/bwa.nix;
      fastqc = callBionix ./tools/fastqc.nix;
      mosdepth = callBionix ./tools/mosdepth.nix;
      platypus = callBionix ./tools/platypus.nix;
      ref = callBionix ./references.nix;
      samtools = callBionix ./tools/samtools.nix;
      strelka = callBionix ./tools/strelka.nix;

      qsub = nixpkgs.callPackage ./lib/qsub.nix {};
  });
in bionix
