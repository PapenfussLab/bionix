{ nixpkgs ? import <nixpkgs> { }, overlays ? [ ] }:

let
  inherit (nixpkgs) fetchurl callPackage;

  bionix = nixpkgs.lib.makeExtensible (self:
    let callBionix = file: attrs: import file ({ bionix = self; } // attrs);
    in
    with self; {
      inherit callBionix;
      id = x: x;
      exec = id;
      exec' = f: exec (_: f) { };
      exec'' = f: exec' (_: f) { };
      callBionixE = p: exec (callBionix p);

      types = callBionix ./lib/types.nix { };
      fetchgdrive = callBionix ./lib/google.nix { };

      bowtie = callBionix ./tools/bowtie.nix { };
      bwa = callBionix ./tools/bwa.nix { };
      cnvkit = callBionix ./tools/cnvkit.nix { };
      compression = callBionix ./tools/compression.nix { };
      crumble = callBionix ./tools/crumble.nix { };
      facets = callBionix ./tools/facets.nix { };
      fastqc = callBionix ./tools/fastqc.nix { };
      gridss = callBionix ./tools/gridss.nix { };
      infercnv = callBionix ./tools/infercnv.nix { };
      kallisto = callBionix ./tools/kallisto.nix { };
      mosdepth = callBionix ./tools/mosdepth.nix { };
      mutect = callBionix ./tools/mutect.nix { };
      minimap2 = callBionix ./tools/minimap2.nix { };
      picard = callBionix ./tools/picard.nix { };
      ref = callBionix ./lib/references.nix { };
      samtools = callBionix ./tools/samtools.nix { };
      sambamba = callBionix ./tools/sambamba.nix { };
      snap = callBionix ./tools/snap.nix { };
      snpeff = callBionix ./tools/snpeff.nix { };
      strelka = callBionix ./tools/strelka.nix { };
      ascat = callBionix ./tools/ascat.nix { };
      fastp = callBionix ./tools/fastp.nix { };
      octopus = callBionix ./tools/octopus.nix { };
      snver = callBionix ./tools/snver.nix { };
      hisat2 = callBionix ./tools/hisat2.nix { };
      xenomapper = callBionix ./tools/xenomapper.nix { };
      manta = callBionix ./tools/manta.nix { };
      delly = callBionix ./tools/delly.nix { };
      lumpy = callBionix ./tools/lumpy.nix { };
      lastal = callBionix ./tools/last.nix { };
      whisper = callBionix ./tools/whisper.nix { };
      star = callBionix ./tools/star.nix { };
      genmap = callBionix ./tools/genmap.nix { };
      subread = callBionix ./tools/subread.nix { };
      hatchet = callBionix ./tools/hatchet.nix { };
      pizzly = callBionix ./tools/pizzly.nix { };
      quip = callBionix ./tools/quip.nix { };
      ampliconarchitect = callBionix ./tools/aa.nix { };

      slurm-run = callPackage ./lib/slurm.nix { };
      slurm-exec = f: x: y:
        slurm-run x (f
          (builtins.removeAttrs x [
            "ppn"
            "mem"
            "walltime"
            "partition"
            "slurmFlags"
            "salloc"
            "srun"
          ])
          y);
      slurm = bionix.extend (self: super: { exec = super.slurm-run; });
      qsub = attrs:
        bionix.extend (self: super:
          with self; rec {
            qsubDefs = {
              ppn = 1;
              mem = 1;
              walltime = "24:00:00";
              tmpDir = "/tmp";
              sleepTime = 60;
              queue = null;
              qsubFlags = null;
              qsubPath = "/usr/bin";
            } // attrs;
            qsub = attrs: (callPackage ./lib/qsub.nix { }) (qsubDefs // attrs);
            exec = f: x: y:
              qsub (builtins.intersectAttrs qsubDefs x) (super.exec f
                (builtins.removeAttrs x (builtins.attrNames qsubDefs))
                y);
          });
      def = f: defs: attrs: f (defs // attrs);

      linkOutputs = x:
        let
          cmds =
            let
              recurse = x:
                if x ? type && x.type == "derivation" then
                  x
                else if builtins.typeOf x == "set" then
                  linkOutputs x
                else
                  abort "linkOutputs: unsupported type";
              link = dst: src: lib.optionalString (src != null) ''
                ln -s ${recurse src} $out/${lib.escapeShellArg dst}
              '';
            in
            ''
              mkdir $out
              ${lib.concatStringsSep "\n" (lib.mapAttrsToList link x)}
            '';
        in
        pkgs.stdenvNoCC.mkDerivation {
          name = "link-outputs";
          nativeBuildInputs = [ pkgs.perl ];
          buildCommand = "exec sh ${pkgs.writeScript "make-links" cmds}";
          passthru.linkInputs = x;
        };

      # Fetching files of specific type
      fetchFastQ = attrs:
        with types;
        tagFiletype (filetype.fq { }) (fetchurl attrs);
      fetchFastA = attrs:
        with types;
        tagFiletype (filetype.fa { }) (fetchurl attrs);
      fetchFastQGZ = attrs:
        with types;
        tagFiletype (filetype.gz (filetype.fq { })) (fetchurl attrs);
      fetchFastAGZ = attrs:
        with types;
        tagFiletype (filetype.gz (filetype.fa { })) (fetchurl attrs);

      # Turn a multi-output derivation into a list of derivations
      outputDrvs = drv: map (o: lib.getAttr o drv) drv.outputs;

      # Export nixpkgs and standard library lib
      pkgs = nixpkgs;
      lib = nixpkgs.lib // {
        inherit types;
        shard = callBionix ./lib/shard.nix { };
        concatMapAttrsStringsSep = s: f: a: with nixpkgs.lib; concatStringsSep s (mapAttrsToList f a);
      };
      stage = x@{ name, stripStorePaths ? true, multicore ? false, ... }:
        (if stripStorePaths then strip else x: x)
          (nixpkgs.stdenvNoCC.mkDerivation (x // {
            name = "bionix-" + name;
            inherit multicore;
          }));
      strip = drv:
        let
          strip-store-paths = nixpkgs.callPackage ./strip-store-paths { };
          stripCommand = ''

            function rewriteOutput {
              find ''${!1} -type f -print0 | xargs -0 -n1 strip-store-paths
            }
            for o in $outputs ; do
              rewriteOutput $o
            done
          '';
        in
        drv.overrideAttrs (attrs:
          { nativeBuildInputs = attrs.nativeBuildInputs or [ ] ++ [ strip-store-paths ]; } // (
            if attrs ? buildCommand then {
              buildCommand = attrs.buildCommand + stripCommand;
            } else {
              fixupPhase = (if attrs ? fixupPhase then attrs.fixupPhase else "")
              + stripCommand;
            }
          ));

      # splitting/joining
      splitFile = file: drv:
        stage {
          name = "split-${file}";
          buildCommand = "ln -s ${drv}/${file} $out";
        };
      split = drv: lib.mapAttrs (p: _: splitFile p drv) (builtins.readDir drv);
      join = drvs:
        stage {
          name = "join";
          buildCommand = ''
            mkdir $out
            ${builtins.concatStringsSep "\n" (builtins.attrValues
              (lib.mapAttrs (n: d: "ln -s ${d} $out/${n}") drvs))}
          '';
        };
      each = f: drv: join (lib.mapAttrs (_: f) (split drv));
    });

  overlayByType = {
    lambda = bionix: overlay:
      bionix.extend
        (self: super: nixpkgs.lib.recursiveUpdate super (overlay self super));
    path = bionix: path: overlay bionix (import path);
  };
  overlay = bionix: overlay:
    let overlayType = builtins.typeOf overlay;
    in
    if overlayByType ? "${overlayType}" then
      overlayByType."${overlayType}" bionix overlay
    else
      builtins.throw ("cannot overlay type " + overlayType);
in
with nixpkgs.lib; foldl overlay bionix overlays
