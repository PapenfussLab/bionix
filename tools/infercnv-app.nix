{stdenv
  ,fetchurl
  ,fetchFromGitHub
  ,callPackage
  ,nettools
  ,rPackages
  ,rWrapper
  ,R
  ,jre
  ,python
  ,pythonPackages
  ,darwin
  ,gettext
  ,gfortran
,makeWrapper}:

let
  buildRPackage = callPackage "${<nixpkgs>}/pkgs/development/r-modules/generic-builder.nix" {
    inherit R gettext gfortran;
    inherit (darwin.apple_sdk.frameworks) Cocoa Foundation;
  };

  GMD = buildRPackage rec {
    name = "GMD-${version}";
    version = "0.3.3";
    src = fetchurl {
      url = "https://cran.r-project.org/src/contrib/Archive/GMD/GMD_0.3.3.tar.gz";
      sha256 = "0jshdcmqcr7lz4p5xb76qbaqavm2609r01lhi9hd0aqnnry18kmg";
    };
    buildInputs = with rPackages; [ R gplots ];
  };

  NGCHMR = buildRPackage rec {
    name = "NGCHMR-${version}";
    version = "git";
    src = fetchFromGitHub {
      owner = "bmbroom";
      repo = "NGCHMR";
      rev = "9f5f1fbf39339d21295b5056e469edcdcbaae142";
      sha256 = "0paw3fz22kbk4ps4mfxzfchqvipspl7a60jsz46fsg10v6d3z7yv";
    };
    propagatedBuildInputs = with rPackages; [ R tsvio digest httr jsonlite nettools ];
    patches = [ ./infercnv-ngchmr.patch ];
  };

  tsvio = buildRPackage rec {
    name = "tsvio-${version}";
    version = "git";
    src = fetchFromGitHub {
      owner = "bmbroom";
      repo = "tsvio";
      rev = "067b01ffc1491d50fc1e104b1fe36208a3997980";
      sha256 = "05byfn2bim51wswffs9lm23p4i0bghyn63rny480dvagydn1a85c";
    };
  };

  inferCNV = buildRPackage rec {
    name = "inferCNV-${version}";
    version = "git";
    requireX = false;
    src = fetchFromGitHub {
      owner = "broadinstitute";
      repo = "inferCNV";
      rev = "cf442af0db6191fa8ba57c4921ac2d1f98c2c39d";
        sha256 = "0cv8qiaqpd6b4152dplnzvgv77cmk961rmvzr27qgmlaazc5hblh";
    };
    propagatedBuildInputs = with rPackages; [ R GMD NGCHMR RColorBrewer gplots optparse logging ];
  };

  r = rWrapper.override {
    packages = with rPackages; [ inferCNV ape ];
  };

  py = python.withPackages (pkgs: with pkgs; [ statistics ]);

  shaidymapgen = fetchurl {
    url = "http://tcga.ngchm.net/NGCHM/ShaidyMapGen.jar";
    sha256 = "1pz710ig8nnydz329ry8fydccbrp3arp614dgba3bcyy9flm3gnw";
  };

in stdenv.mkDerivation rec {
  name = inferCNV.name;
  src = inferCNV.src;
  buildInputs = [ r makeWrapper py ];
  propagatedBuildInputs = [ jre ];
  installPhase = ''
    mkdir -p $out/bin
    cp scripts/* $out/bin
  '';
  postFixup = ''
    wrapProgram $out/bin/inferCNV.R --set SHAIDYMAPGEN=${shaidymapgen}
  '';
}
