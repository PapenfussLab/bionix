{ buildRPackage
, fetchFromGitHub
, R
, htslib
, zlib
}:

let
  pctGCdata = buildRPackage rec {
    name = "pctGCdata-${version}";
    version = "0.3.0";
    requireX = false;
    src = fetchFromGitHub {
      owner = "mskcc";
      repo = "pctGCdata";
      rev = "v0.3.0";
      sha256 = "sha256-g0/TlV/a23UuogBc/0m6QTDaHXxFGtmdOSQTiPRGNbE=";
    };
    buildInputs = [ R ];
  };

in
buildRPackage rec{
  name = "facets-${version}";
  version = "0.6.1";
  requireX = false;
  src = fetchFromGitHub {
    owner = "mskcc";
    repo = "facets";
    rev = "v${version}";
    sha256 = "sha256-GbEvhdkS+lCcJclCNtdXJgku6+WhXsUiSp4p2uhtqUs=";
  };
  buildInputs = [ R htslib zlib ];
  propagatedBuildInputs = [ pctGCdata ];
  postBuild = ''
    cd inst/extcode
    g++ --std=c++11 snp-pileup.cpp -lhts -o snp-pileup
    cd ../../
  '';
  postInstall = ''
    mkdir -p $out/bin
    cp inst/extcode/snp-pileup $out/bin
  '';
}
