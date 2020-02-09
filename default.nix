{nixpkgs ? import <nixpkgs> {}}:

let
  inherit (nixpkgs) fetchurl callPackage;

  bionix = nixpkgs.lib.makeExtensible (self:
  let callBionix = file: attrs: import file ({ bionix = self; } // attrs);
  in with self; {
    callBionix = callBionix;
    id = x: x;
    exec = f: x: y: f x y;
    exec' = f: exec (_: f) {};
    exec'' = f: exec' (_: f) {};
    callBionixE = p: exec (callBionix p);

    types = callBionix ./lib/types.nix {};

    bowtie = callBionix ./tools/bowtie.nix {};
    bwa = callBionix ./tools/bwa.nix {};
    cnvkit = callBionix ./tools/cnvkit.nix {};
    compression = callBionix ./tools/compression.nix {};
    crumble = callBionix ./tools/crumble.nix {};
    facets = callBionix ./tools/facets.nix {};
    fastqc = callBionix ./tools/fastqc.nix {};
    gridss = callBionix ./tools/gridss.nix {};
    infercnv = callBionix ./tools/infercnv.nix {};
    kallisto = callBionix ./tools/kallisto.nix {};
    mosdepth = callBionix ./tools/mosdepth.nix {};
    mutect = callBionix ./tools/mutect.nix {};
    minimap2 = callBionix ./tools/minimap2.nix {};
    picard = callBionix ./tools/picard.nix {};
    platypus = callBionix ./tools/platypus.nix {};
    ref = callBionix ./lib/references.nix {};
    samtools = callBionix ./tools/samtools.nix {};
    snap = callBionix ./tools/snap.nix {};
    snpeff = callBionix ./tools/snpeff.nix {};
    strelka = callBionix ./tools/strelka.nix {};
    ascat = callBionix ./tools/ascat.nix {};
    fastp = callBionix ./tools/fastp.nix {};
    octopus = callBionix ./tools/octopus.nix {};
    snver = callBionix ./tools/snver.nix {};
    hisat2 = callBionix ./tools/hisat2.nix {};
    xenomapper = callBionix ./tools/xenomapper.nix {};
    manta = callBionix ./tools/manta.nix {};
    delly = callBionix ./tools/delly.nix {};
    lumpy = callBionix ./tools/lumpy.nix {};
    lastal = callBionix ./tools/last.nix {};

    slurm = attrs: bionix.extend (self: super: with self; rec {
      slurmDefs = { ppn = 1; mem = 1; walltime = "24:00:00"; partition = null; slurmFlags = null; salloc = "/usr/bin/salloc"; srun = "/usr/bin/srun"; } // attrs;
      slurm = attrs: (callPackage ./lib/slurm.nix {}) (slurmDefs // attrs);
      exec = f: x: y: slurm (builtins.intersectAttrs slurmDefs x) (super.exec f (builtins.removeAttrs x (builtins.attrNames slurmDefs)) y);
      });
    qsub = attrs: bionix.extend (self: super: with self; rec {
      qsubDefs = { ppn = 1; mem = 1; walltime = "24:00:00"; tmpDir = "/tmp"; sleepTime = 60; queue = null; qsubFlags = null;  qsubPath = "/usr/bin"; } // attrs;
      qsub = attrs: (callPackage ./lib/qsub.nix {}) (qsubDefs // attrs);
      exec = f: x: y: qsub (builtins.intersectAttrs qsubDefs x) (super.exec f (builtins.removeAttrs x (builtins.attrNames qsubDefs)) y);
    });
    def = f: defs: attrs: f (defs // attrs);
    pipe = let g = fs: with builtins; let h = head fs; t = tail fs; in if t != [] then x: (g t (h x)) else h; in g;

    linkOutputs = x: with lib; nixpkgs.stdenvNoCC.mkDerivation {
      name = "link-outputs";
      outputs = [ "out" ] ++ attrNames x;
      nativeBuildInputs = [ pkgs.perl ];
      buildCommand = let
        recurse = x: if x ? type && x.type == "derivation" then x else
          if builtins.typeOf x == "set" then linkOutputs x
          else error "linkOutputs: unsupported type";
        link = dst: src: ''
          ln -s ${recurse src} $(perl -e 'print $ENV{"${dst}"}') ; ln -s ${recurse src} $out/${dst}
        '';
      in "mkdir $out \n" + (concatStringsSep "\n" (mapAttrsToList link x));
    };

    # Fetching files of specific type
    fetchFastQ = attrs: with types; tagFiletype (filetype.fq {}) (fetchurl attrs);
    fetchFastA = attrs: with types; tagFiletype (filetype.fa {}) (fetchurl attrs);
    fetchFastQGZ = attrs: with types; tagFiletype (filetype.gz (filetype.fq {})) (fetchurl attrs);
    fetchFastAGZ = attrs: with types; tagFiletype (filetype.gz (filetype.fa {})) (fetchurl attrs);

    # Turn a multi-output derivation into a list of derivations
    outputDrvs = drv: map (o: lib.getAttr o drv) drv.outputs;

    # Export nixpkgs and standard library lib
    pkgs = nixpkgs;
    lib = nixpkgs.lib // { types = types; shard = callBionix ./lib/shard.nix {};};
    stage = x@{ name, stripStorePaths ? true, multicore ? false, ... }:
      (if stripStorePaths then strip else x: x) (nixpkgs.stdenvNoCC.mkDerivation (x // {name = "bionix-" + name; inherit multicore;}));
    strip = drv: let
        stripCommand = ''

        function rewrite {
          sed -i 's|[A-Za-z0-9+/]\{32\}-bionix|00000000000000000000000000000000-bionix|g' $1
        }
        function rewriteOutput {
          if [ -f ''${!1} ] ; then
            rewrite ''${!1}
          else
            for f in $(find ''${!1} -type f) ; do
              rewrite $f
            done
          fi
        }
        for o in $outputs ; do
          rewriteOutput $o
        done
      '';
      in drv.overrideAttrs (attrs:
        if attrs ? buildCommand then
          {buildCommand = attrs.buildCommand + stripCommand;}
        else
          { fixupPhase = (if attrs ? fixupPhase then attrs.fixupPhase else "") + stripCommand; }
        );

    # splitting/joining
    splitFile = file: drv: stage {
      name = "split-${file}";
      buildCommand = "ln -s ${drv}/${file} $out";
    };
    split = drv: lib.mapAttrs (p: _: splitFile p drv) (builtins.readDir drv);
    join = drvs: stage {
      name = "join";
      buildCommand = ''
        mkdir $out
        ${builtins.concatStringsSep "\n" (builtins.attrValues (lib.mapAttrs (n: d: "ln -s ${d} $out/${n}") drvs))}
      '';
    };
    each = f: drv: join (lib.mapAttrs (_: f) (split drv));
  });
in bionix
