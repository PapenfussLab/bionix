{ python3Packages, fetchFromGitHub, cmake, gurobi, writeText }:

let

  findgurobi = writeText "FindGUROBI.cmake" ''
    set(GUROBI_CPP_LIB ${gurobi}/lib/libgurobi_c++.a)
    set(GUROBI_LIB ${gurobi}/lib/libgurobi91.so)
    set(GUROBI_INCLUDE_DIR ${gurobi}/include)
    set(GUROBI_LIBRARIES ''${GUROBI_CPP_LIB} ''${GUROBI_LIB} -lpthread)
    set(GUROBI_FOUND TRUE)
  '';

in
python3Packages.buildPythonApplication rec {
  pname = "HATCHet";
  version = "0.4.9";

  src = fetchFromGitHub {
    owner = "raphael-group";
    repo = "hatchet";
    rev = "v${version}";
    sha256 = "sha256-MB9XFbkLQTf6ZUPrisSzGU8Jeq6SrlMMCQtoyvx/Xvc=";
  };

  dontConfigure = true;

  patchPhase = ''
    cat ${findgurobi} > FindGUROBI.cmake
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
  ];
}
