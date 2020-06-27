{ stdenv, fetchFromGitHub, fetchurl }:

stdenv.mkDerivation {
  pname = "whisper";
  version = "2.0";

  src = fetchFromGitHub {
    owner = "refresh-bio";
    repo = "Whisper";
    rev = "v2.0";
    sha256 = "100p2pqhli123wnqkxrvjcwnmlvcp1rk63whh121jannnw81rh2m";
  };

  patches = [ (fetchurl {
    url = "https://github.com/refresh-bio/Whisper/commit/9ad77e7de68d91d9427cfbe6211e83cea89206ab.patch";
    sha256 = "0i1s70rq5z4b3yik76vs1cxri4n3mpgfn91k35b3mwzhpxvfh7nn";}) ];

  preConfigure = ''
    cd src

    # disable default static linking
    sed -i 's/ -static / /' makefile
  '';

  installPhase = ''
    mkdir -p $out/bin
    cp whisper whisper-index $out/bin
  '';

  buildPhase = "make -j $NIX_BUILD_CORES";

  meta = with stdenv.lib; {
    description = "Short read sequence mapper";
    license = licenses.gpl3;
    homepage = "https://github.com/refresh-bio/whisper";
    maintainers = with maintainers; [ jbedo ];
    platforms = platforms.linux;
  };

}

