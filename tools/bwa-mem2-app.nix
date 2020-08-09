{ clangStdenv, fetchFromGitHub, zlib }:
  clangStdenv.mkDerivation rec {
  pname = "bwa-mem2";
  version = "2.0";
  src = fetchFromGitHub {
    fetchSubmodules = true;
    owner = "bwa-mem2";
    repo = "bwa-mem2";
    rev = "v${version}";
    sha256 = "0q5wqal0nfxd3yfbmxahyaiqqmsrrplnwhplcjvz1xzw7bxwwnnj";
  };
  buildInputs = [ zlib ];
  installPhase = ''
    mkdir -p $out/bin
    cp bwa-mem2 $out/bin
  '';
  disableHardening = [ "pie" ];
}
