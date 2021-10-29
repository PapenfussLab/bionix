{ bionix }:

{ url, sha256 }:

with bionix.pkgs;

let
  gdown = fetchurl {
    url =
      "https://raw.githubusercontent.com/circulosmeos/gdown.pl/master/gdown.pl";
    sha256 = "1pw3vg70bgf33akbbphpr6zn3jndv0khmsa3k0m877hgzg1v52qv";
  };
in
runCommand "gdown"
{
  nativeBuildInputs = [ perl wget ];
  outputHashAlgo = "sha256";
  outputHash = sha256;
  outputHashMode = "flat";
} ''
  perl ${gdown} "${url}" $out
''
