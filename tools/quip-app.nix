{
  stdenv,
  fetchFromGitHub,
  autoreconfHook,
  zlib,
}:
stdenv.mkDerivation rec {
  name = "quip";

  src = fetchFromGitHub {
    owner = "dcjones";
    repo = "quip";
    rev = "9165bb5ac17a2cf2822e437df573487d7adfa1cc";
    sha256 = "PpvgLzdiYcLyJ1JH9UzbUd29eKIlyY440HiJGYWvkrs=";
  };

  nativeBuildInputs = [autoreconfHook];
  buildInputs = [zlib];
}
