{ clangStdenv, fetchFromGitHub, zlib }:
clangStdenv.mkDerivation rec {
  pname = "bwa-mem2";
  version = "2.2.1";

  src = fetchFromGitHub {
    fetchSubmodules = true;
    owner = "bwa-mem2";
    repo = "bwa-mem2";
    rev = "v${version}";
    sha256 = "sha256-2wDhtTlxIsJAlQ+Z72kmwbVqD9GZnd845eIRYMxOzgU=";
  };

  buildInputs = [ zlib ];

  installPhase = ''
    mkdir -p $out/bin
    cp bwa-mem2 bwa-mem2.* $out/bin
  '';

  NIX_CFLAGS_COMPILE = "-Wno-register -Wno-implicit-function-declaration";
}
