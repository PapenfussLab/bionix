{ stdenv
, mosdepth
, python
, flags ? null}:

input:

stdenv.mkDerivation {
  name = "mosdepth-depth";
  buildInputs = [ python ];
  buildCommand = ''
    python ${mosdepth.src}/scripts/plot-dist.py ${input}/*global.dist.txt
    mv dist.html $out
  '';
}
