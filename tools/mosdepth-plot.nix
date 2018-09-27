{ stdenv
, lib
, mosdepth
, python
, flags ? null}:

with lib;

{ inputs
, names ? []}:

stdenv.mkDerivation {
  name = "mosdepth-depth";
  buildInputs = [ python ];
  buildCommand = ''
    python ${mosdepth.src}/scripts/plot-dist.py ${concatMapStringsSep " " (i: "${i}/out.mosdepth.global.dist.txt") inputs}
    ${concatStringsSep "\n" (zipListsWith (i: n: "substituteInPlace dist.html --replace ${i}/out ${n}") inputs names)}
    substituteInPlace dist.html --replace "x: 0.1" "x: 0.9"
    substituteInPlace dist.html --replace "y: 0.1" "y: 0.9"
    mv dist.html $out
  '';
}
