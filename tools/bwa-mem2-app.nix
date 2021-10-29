{ clangStdenv, fetchFromGitHub, zlib }:
clangStdenv.mkDerivation rec {
  pname = "bwa-mem2";
  version = "2.1";
  src = fetchFromGitHub {
    fetchSubmodules = true;
    owner = "bwa-mem2";
    repo = "bwa-mem2";
    rev = "v${version}";
    sha256 = "sha256-T0nkO+NehAMFuwGi7HrQxq57gEsA16TlEWa2SHerxV4=";
  };
  buildInputs = [ zlib ];
  installPhase = ''
    mkdir -p $out/bin
    cp bwa-mem2 bwa-mem2.* $out/bin
  '';
  disableHardening = [ "pie" ];
}
