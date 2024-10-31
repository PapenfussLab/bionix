{ stdenv, zig }:

stdenv.mkDerivation {
  name = "strip-store-paths";
  nativeBuildInputs = [ zig ];
  src = ./strip-store-paths.zig;
  
  XDG_CACHE_HOME = "Cache";

  unpackPhase = ''
    cp $src strip-store-paths.zig
  '';

  buildPhase = ''
    zig build-exe -OReleaseFast -mcpu=baseline strip-store-paths.zig
  '';

  installPhase = ''
    install -Dm 755 ./strip-store-paths $out/bin/strip-store-paths
  '';
}
