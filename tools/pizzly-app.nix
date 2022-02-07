{ lib, stdenv, fetchFromGitHub, cmake, zlib }:

stdenv.mkDerivation rec {
  pname = "pizzly";
  version = "0.37.3";

  src = fetchFromGitHub {
    repo = pname;
    owner = "pmelsted";
    rev = "v${version}";
    sha256 = "sha256-VB3kYjnhg9+QMdCAbYAzf30kO8zXZlEvK913ALZiUwU=";
  };

  nativeBuildInputs = [ cmake ];
  buildInputs = [ zlib ];

  installPhase = ''
    runHook preInstall
    install -Dm 755 pizzly $out/bin/pizzly
    runHook postInstall
  '';

  meta = with lib; {
    description = "Fusion detection from kallisto";
    homepage = "https://github.com/pmelsted/pizzly";
    license = licenses.bsd2;
    platforms = platforms.linux;
    maintainers = with maintainers; [ jbedo ];
  };
}
