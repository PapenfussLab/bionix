{ lib
, fetchurl
, rPackages
, rWrapper
, buildPythonPackage
, fetchPypi
, biopython
, six
, numpy
, scipy
, pandas
, matplotlib
, reportlab
, zlib
, bzip2
, htslib
, libjpeg
, pkgconfig
, futures
}:

let pyfaidx = buildPythonPackage rec {
      pname = "pyfaidx";
      version = "0.5.3.1";

      src = fetchPypi {
        inherit pname version;
        sha256 = "0mjbksbj9hh2cf0yjr951cjahhn0lg7p71kd3kvbnscqyxa44kfr";
      };

      propagatedBuildInputs = [ six ];
    };

    pysam = buildPythonPackage rec {
      pname = "pysam";
      version = "0.14";

      src = fetchPypi {
        inherit pname version;
        sha256 = "11snrwl5pn6knarf0c7zcpkhc7rijb6cfpw8sl3pdvnynj1fmj69";
      };

      buildInputs = [ zlib ];

      preConfigure = ''
        export HTSLIB_MODE="external"
        export HTSLIB_LIBRARY_DIR=${htslib}/lib
        export HTSLIB_INCLUDE_DIR=${htslib}/include
      '';

    };

    future = buildPythonPackage rec {
      pname = "future";
      version = "0.16.0";

      doCheck = false;

      src = fetchPypi {
        inherit pname version;
        sha256 = "1nzy1k4m9966sikp0qka7lirh8sqrsyainyf8rk97db7nwdfv773";
      };
    };

    pillow = buildPythonPackage rec {
      pname = "Pillow";
      version = "5.0.0";

      doCheck = false;

      src = fetchPypi {
        inherit pname version;
        sha256 = "1fz1n1cq65dqdbiwy1cn8l21amqbzq18rdmmcr670ks24dn9vwhj";
      };

      buildInputs = [ zlib libjpeg pkgconfig ];
    };

    cghFLasso = rPackages.buildRPackage rec {
      name = "cghFLasso-${version}";
      version = "0.2-1";
      src = fetchurl {
        url = "https://cran.r-project.org/src/contrib/Archive/cghFLasso/cghFLasso_${version}.tar.gz";
        sha256 = "0b1hnjf9g0v47hbz0dy9m6jhcl1ky20yyhhmm8myng2sndcpjsbf";
      };
    };

    cnvR = rWrapper.override {
      packages = with rPackages; [ DNAcopy cghFLasso ];
    };

in buildPythonPackage rec {
  pname = "CNVkit";
  version = "0.9.5";

  src = fetchPypi {
    inherit pname version;
    sha256 = "1sa70bmnxj1lzp33pbj3axk6n77czswwj9cirimxh2qrn84i7vs3";
  };

  propagatedBuildInputs = [ biopython numpy scipy pandas matplotlib reportlab pyfaidx pysam futures future pillow cnvR ];
}
