{stdenv, fetchFromGitHub, autoreconfHook, htslib}:

stdenv.mkDerivation rec {
  name = "crumble-${version}";
  version = "git";

  src = fetchFromGitHub {
    owner = "jkbonfield";
    repo = "crumble";
    rev = "217fa927a5fffe7dfeed9ebf25ae78d5b513afac";
    sha256 = "1sjgbz0dj9pszxqgkrslcdx305s50a5gyx3f5ajsnm4bynsjmv9i";
  };

  buildInputs = [ autoreconfHook htslib ];
}
