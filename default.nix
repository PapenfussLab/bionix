{nixpkgs ? import <nixpkgs> {}}:

let
  bionix = nixpkgs.lib.makeExtensible (self:
  let callBionix = file: attrs: import file ({ bionix = self; nixpkgs = nixpkgs; } // attrs);
  in with self; {
    callBionix = callBionix;
    id = x: x;

    types = callBionix ./lib/types.nix {};

    bwa = callBionix ./tools/bwa.nix {};
    compression = callBionix ./tools/compression.nix {};
    crumble = callBionix ./tools/crumble.nix {};
    fastqc = callBionix ./tools/fastqc.nix {};
    gridss = callBionix ./tools/gridss.nix {};
    infercnv = callBionix ./tools/infercnv.nix {};
    mosdepth = callBionix ./tools/mosdepth.nix {};
    mutect = callBionix ./tools/mutect.nix {};
    platypus = callBionix ./tools/platypus.nix {};
    samtools = callBionix ./tools/samtools.nix {};
    strelka = callBionix ./tools/strelka.nix {};

    qsub = nixpkgs.callPackage ./lib/qsub.nix {};
    qsubAttr = qsubAttrs: f: attrs: i: qsub qsubAttrs (f attrs i);
    qsubAttrs = attrs: nixpkgs.lib.mapAttrs (_: x: qsubAttr attrs x);
    ref = callBionix ./lib/references.nix {};
    def = f: defs: attrs: f (defs // attrs);
    defQsub = qsubAttrs: f: defs: qsubAttr qsubAttrs (def f defs);

  });
in bionix
