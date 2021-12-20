{ stdenv, zig }:

stdenv.mkDerivation {
  name = "strip-store-paths";
  nativeBuildInputs = [ zig ];
  src = ./strip-store-paths.zig;

  unpackPhase = ''
    cp $src strip-store-paths.zig
  '';

  buildPhase = ''
    export HOME=$TMPDIR
    zig build-exe -OReleaseFast strip-store-paths.zig
  '';

  installPhase = ''
    install -Dm 755 ./strip-store-paths $out/bin/strip-store-paths
  '';
}
