{stdenv, fetchurl, makeWrapper, unzip, fetchFromGitHub}:

let
  oldnix = import (fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs";
    rev = "83a893c38a83877588e3ca7ccfeabaa973c30acd";
    sha256 = "0q7214hag7h95irvhkdb648m09b9jspb0raw1qjrx7y4grzb165h";
  }) {};

  jre = oldnix.openjdk7;

in stdenv.mkDerivation rec {
  name = "mutect-${version}";
  version = "1.1.5";

  src = fetchurl {
    url = "https://github.com/broadinstitute/mutect/releases/download/${version}/muTect-${version}-bin.zip";
    sha256 = "1pq7iv720bp970qsyyshwk98xdb7naw566y6gk9cpj6bmm08z9v3";
  };

  buildInputs = [ makeWrapper unzip ];

  unpackPhase = ''
    unzip $src -d $TMPDIR
  '';
  installPhase = ''
    install -Dt $out/libexec/mutect muTect-${version}.jar
    mkdir $out/bin
    makeWrapper ${jre}/bin/java $out/bin/mutect --add-flags "-jar $out/libexec/mutect/muTect-${version}.jar"
  '';
}
