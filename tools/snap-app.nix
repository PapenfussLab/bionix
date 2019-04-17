{stdenv, fetchFromGitHub, zlib}:

stdenv.mkDerivation rec {
  name = "snap-git";
  src = fetchFromGitHub {
    owner = "amplab";
    repo = "snap";
    rev = "669a341cfa5c3e0a80098bf5aa51e3332f868d2b";
    sha256 = "1810d8wn88wbagip7hlnwwg7z6jqyiyfnnvdiir06s3rrmn95mjh";
  };
  buildInputs = [ zlib ];
  postConfigure = ''
    sed -i 's/-Wno-format/-Wformat/g' Makefile
  '';
  installPhase = ''
    mkdir -p $out/bin
    cp snap-aligner $out/bin
  '';
}
