{
  stdenv,
  fetchurl,
  fetchFromGitHub,
  python3,
}: let
  python = python3.withPackages (pkgs:
    with pkgs; [
      numpy
      scipy
      pysam
      matplotlib
      future
      (mosek pkgs)
    ]);

  mosek = assert stdenv.system == "x86_64-linux";
    pkgs:
      pkgs.buildPythonPackage {
        pname = "mosek";
        version = "8.1.0.83";
        src = fetchurl {
          url = "https://download.mosek.com/stable/8.1.0.83/mosektoolslinux64x86.tar.bz2";
          sha256 = "sha256-d/S/IalmQwizWYZ89ZskUoVAaXWYszuw7w+w0Vp+13k";
        };
        doCheck = false;
        preBuild = ''
          cd 8/tools/platform/linux64x86/python/3/
        '';
        propagatedBuildInputs = with pkgs; [numpy];
        postInstall = ''
          find $out -name lib\*.so\* -print0 | xargs -0 \
            patchelf --add-rpath ${stdenv.cc.cc.lib}/lib
        '';
      };
in
  stdenv.mkDerivation rec {
    pname = "AmpliconArchitect";
    version = "1.3";

    src = fetchFromGitHub {
      owner = "virajbdeshpande";
      repo = pname;
      rev = "40da8520a953810ad43e5a6fdf4aba7449d7f5e0";
      sha256 = "sha256-4SAOpdjXiZFTfpD6WpLfs2zDyGT2hcWabl+sUjboBpc=";
    };

    doBuild = false;
    installPhase = ''
      mkdir -p $out/libexec
      cp -r src $out/libexec/aa
      mkdir $out/bin
      ln -s $out/libexec/aa/{AmpliconArchitect,amplified_intervals,ref_util,downsample}.py $out/bin
    '';

    buildInputs = [python];
  }
