{ stdenv, unzip, perl, fetchurl, jre }:

stdenv.mkDerivation rec {
  name = "fastqc-${version}";
  version = "0.11.7";

  src = fetchurl {
    url = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip";
    sha256 = "04hifbfrh60s5kzqr7n46bcviaiymr1rx67b88s3cpxydf3m1ksr";
  };

  buildInputs = [ unzip perl ];

  phases = [ "unpackPhase" "installPhase" "fixupPhase" ];

  installPhase = ''
    mkdir -p $out/libexec/fastqc
    cp -r . $out/libexec/fastqc
    mkdir $out/bin
    ln -s $out/libexec/fastqc/fastqc $out/bin
    substituteInPlace $out/bin/fastqc --replace "my \$java_bin = 'java'" "my \$java_bin = '${jre}/bin/java'"
    chmod 755 $out/libexec/fastqc/fastqc
  '';
}
