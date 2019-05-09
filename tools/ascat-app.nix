{ stdenv
  , callPackage
  , buildPerlPackage
  , fetchurl
  , fetchFromGitHub
  , perlPackages
  , rWrapper
  , rPackages }:
  with callPackage ./bioperl.nix {};
  let
    ascat = fetchurl {
      url = "https://raw.githubusercontent.com/Crick-CancerGenomics/ascat/v2.5.1/ASCAT/R/ascat.R";
      sha256 = "1rja9s6rksmi0kc6lhx1vr5yqv2xazgxdlwc7mbj2m881x8nngb1";
      };
    ascatR = rWrapper.overrideAttrs (attrs:
      {
        packages = with rPackages;
        [ RColorBrewer ];
        });
    ascatNGS = buildPerlPackage rec {
      name = "ascatNGS-${version}";
      version = "4.1.2";
      src = fetchFromGitHub {
        owner = "cancerit";
        repo = "ascatNgs";
        rev = "v${version}";
        sha256 = "0qc1wp37k2c9vya992cfiphp2jpp8hvfvqz5ydi2d5b71prfaw82";
        };
      preConfigure = ''
        cd perl
        cp ${ascat} share/ascat/ascat.R
        '';
      nativeBuildInputs = with perlPackages;
      [ TestFatal TestWarn ];
      propagatedBuildInputs = with perlPackages;
      [
        ascatR
        FileShareDirInstall
        pcapCore
        FileShareDir
        cgpVcf
        alleleCountPerl
        ];
      };
  in ascatNGS
