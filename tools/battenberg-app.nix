{ stdenv
,callPackage
  , fetchFromGitHub
  , perl
  , perlPackages
  , buildPerlPackage
  }:
  with callPackage ./bioperl.nix {};
  buildPerlPackage rec {
    name = "cgpBattenberg-${version}";
    version = "3.3.0";
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "cgpBattenberg";
      rev = "v${version}";
      sha256 = "19xdpnsgx2yw99lgyrkw3pq943j2kid5xb0qrfzs9wax98303pxg";
      };
    nativeBuildInputs = with perlPackages;
    [ TryTiny ];
    propagatedBuildInputs = with perlPackages;
    [
      FileShareDirInstall
      ConstFast
      IPCSystemSimple
      FileWhich
      cgpVcf
      ListMoreUtils
      pcapCore
      ];
    preConfigure = "cd perl";
    doCheck = false;
    }
