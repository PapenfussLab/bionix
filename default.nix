{nixpkgs ? import <nixpkgs> {}}:

let
  bionix = nixpkgs.lib.makeExtensible (self:
  let callBionix = file: attrs: import file ({ bionix = self; nixpkgs = nixpkgs; } // attrs);
  in with self; {
      callBionix = callBionix;
      id = x: x;

      bwa = callBionix ./tools/bwa.nix {};
      crumble = callBionix ./tools/crumble.nix {};
      fastqc = callBionix ./tools/fastqc.nix {};
      gridss = callBionix ./tools/gridss.nix {};
      mosdepth = callBionix ./tools/mosdepth.nix {};
      platypus = callBionix ./tools/platypus.nix {};
      samtools = callBionix ./tools/samtools.nix {};
      strelka = callBionix ./tools/strelka.nix {};

      qsub = nixpkgs.callPackage ./lib/qsub.nix {};
      qsubAttr = qsubAttrs: f: attrs: i: qsub qsubAttrs (f attrs i);
      qsubAttrs = attrs: nixpkgs.lib.mapAttrs (_: x: qsubAttr attrs x);
      ref = callBionix ./lib/references.nix {};
  });
in bionix
