{ stdenv
, fetchFromGitHub
, zlib 
}:

stdenv.mkDerivation rec {
  name = "fastp-${version}";
  version = "0.20.0";

  src = fetchFromGitHub {
    owner = "OpenGene";
    repo = "fastp";
    rev = "v${version}";
    sha256 = "0y0qfp3j3gqnmlqskna8x43acss21vxwck287c4fagxlcaba0s30";
  };

  buildInputs = [ zlib ];

  installPhase = ''
    mkdir -p $out/bin
    cp fastp $out/bin
  '';
}
