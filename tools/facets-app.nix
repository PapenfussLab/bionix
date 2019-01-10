{buildRPackage
,fetchFromGitHub
,R
,htslib
,zlib}:

let
  pctGCdata = buildRPackage rec {
    name = "pctGCdata-${version}";
    version = "0.2.0";
    requireX = false;
    src = fetchFromGitHub {
      owner = "mskcc";
      repo = "pctGCdata";
      rev = "v${version}";
      sha256 = "1qq0fmm3zwz6rv0ka82850ww0qj50621gln9i0gfs8k3wyqil4l8";
    };
    buildInputs = [ R ];
  };

in buildRPackage rec{
  name = "facets-${version}";
  version = "0.5.6";
  requireX = false;
  src = fetchFromGitHub {
    owner = "mskcc";
    repo = "facets";
    rev = "v${version}";
    sha256 = "0yqr23446y2h8xgbgj0r6sl4i778111drgsi06a2cqy530xcmbxs";
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
