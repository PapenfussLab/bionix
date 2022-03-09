{ lib, runCommand, python3Packages, fetchFromGitHub, cmake, gurobi, writeText }:

let

  gurobi' = runCommand "gurobi${lib.replaceStrings ["."] [""] gurobi.version}" {} "ln -s ${gurobi} $out";
  in

python3Packages.buildPythonApplication rec {
  pname = "HATCHet";
  version = "0.4.14";

  src = fetchFromGitHub {
    owner = "raphael-group";
    repo = "hatchet";
    rev = "v${version}";
    sha256 = "sha256-zkbjwdtvRNsZWvhtQy8TA3o68l7Uf4WQby1/3/sHq98=";
  };

  dontConfigure = true;
  GUROBI_HOME=gurobi';

  patchPhase = ''
    sed -i 's/''${GUROBI_LIB}/''${GUROBI_LIB} -lpthread/' FindGUROBI.cmake
  '';

  nativeBuildInputs = [ cmake ];
  propagatedBuildInputs = with python3Packages; [
    biopython
    matplotlib
    pandas
    psutil
    pyomo
    pysam
    requests
    seaborn
    scikit-learn
    scipy
    gurobipy
  ];
}
