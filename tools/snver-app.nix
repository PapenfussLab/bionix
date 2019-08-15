{ stdenv, fetchurl, jre, makeWrapper }:

stdenv.mkDerivation rec {
  pname = "SNVer";
  version = "0.5.3";

  src = fetchurl {
    url = "mirror://sourceforge/snver/SNVer-${version}.tar.gz";
    sha256 = "1y3c3gm1zdh4iz6zh1lyaaq1ks205wjm3vwx6wdsnh896xrphf5c";
  };

  buildInputs = [ makeWrapper ];

  unpackPhase = ''
    mkdir -p $out/libexec/SNVer
    tar -zxvf $src -C $out/libexec/SNVer
    rm -rf $out/libexec/SNVer/test
  '';

  installPhase = ''
    mkdir $out/bin
    makeWrapper ${jre}/bin/java $out/bin/SNVerPool --add-flags "-jar $out/libexec/SNVer/SNVerPool.jar"
  '';
}
