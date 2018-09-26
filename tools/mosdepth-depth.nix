{ stdenv
, mosdepth
, flags ? null}:

input:

stdenv.mkDerivation {
  name = "mosdepth-depth";
  buildInputs = [ mosdepth ];
  buildCommand = ''
    mkdir $out
    mosdepth -t $NIX_BUILD_CORES ${optionalString (flags != null) flags} $out/out ${input}
  '';
}
