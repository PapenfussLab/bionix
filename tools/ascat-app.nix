{stdenv, callPackage, buildPerlPackage, fetchurl, fetchFromGitHub, perlPackages, R, bwa, samtools, pkgconfig, zlib, htslib, curl, bzip2, lzma, gnutls, nettle, gmp, p11-kit, libtasn1, perl, psmisc, time, vcftools, rWrapper, rPackages}:

let
  ascat = fetchurl {
    url = "https://raw.githubusercontent.com/Crick-CancerGenomics/ascat/v2.5.1/ASCAT/R/ascat.R";
    sha256 = "1rja9s6rksmi0kc6lhx1vr5yqv2xazgxdlwc7mbj2m881x8nngb1";
  };

  libmaus2 = stdenv.mkDerivation rec {
    name = "libmaus-${version}";
    version = "2.0.499-release-20180606122508";
    src = fetchFromGitHub {
      owner = "gt1";
      repo = "libmaus2";
      rev = version;
      sha256 = "1rk2f3jirn84vgsfml4c8pr6vl5s0xbds1d0givb5pvjmc6jz61x";
    };
    buildInputs = [ zlib ];
  };

  biobambam2 = stdenv.mkDerivation rec {
    name = "biobambam-${version}";
    version = "2.0.87-release-20180301132713";
    src = fetchFromGitHub {
      owner = "gt1";
      repo = "biobambam2";
      rev = version;
      sha256 = "1hjh6kl62hr9afi6fa6x563823lmmmijh1y2204db0nfnaybnrn6";
    };
    nativeBuildInputs = [ pkgconfig ];
    buildInputs = [ libmaus2 zlib ];
  };

  libBigWig = stdenv.mkDerivation rec {
    name = "libBigWig-${version}";
    version = "0.4.2";
    src = fetchFromGitHub {
      owner = "dpryan79";
      repo = "libBigWig";
      rev = version;
      sha256 = "0h2smg24v5srdcqzrmz2g23cmlp4va465mgx8r2z571sfz8pv454";
    };
    buildInputs = [ curl ];
    makeFlags = [ "prefix=$(out)" ];
  };

  cgpBigWig = stdenv.mkDerivation rec {
    name = "cgpBigWig-${version}";
    version = "1.0.2";
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "cgpBigWig";
      rev = version;
      sha256 = "1piprbibmrwbh524plp9s35z5cf1yk7713dxsdik2qv363p6pfsa";
    };
    buildInputs = [ htslib libBigWig curl bzip2 lzma gnutls libtasn1 nettle gmp p11-kit ];
    configurePhase = ''
      cd c
      patchShebangs .
    '';
    makeFlags = [ "HTSLIB=${htslib}" "LIBBIGWIG=${libBigWig}" ];
    installPhase = ''
      mkdir -p $out
      cp -r ../bin $out
    '';
  };

  BioDBHTS = buildPerlPackage rec {
    name = "Bio-DB-HTS-${version}";
    version = "2.11";
    src = fetchFromGitHub {
      owner = "Ensembl";
      repo = "Bio-DB-HTS";
      rev = version;
      sha256 = "0n60vmsd39pz6zbkwz78fdzsmfg8dxgwx0m166vm8515q446922c";
    };
    nativeBuildInputs = with perlPackages; [ ModuleBuild ];
    buildInputs = [ htslib zlib ];
    #propagatedBuildInputs = [ BioRootVersion ];
    configurePhase = ''
      export HTSLIB_DIR=${htslib}
      sed -i 's|"-Wl,-rpath,\$hts_lib",||' Build.PL
      perl ./Build.PL --install_base=$out
    '';
    buildPhase = ''
      ./Build
    '';
    installPhase = ''
      ./Build install
      mkdir -p $out/lib/perl5/site_perl/${perl.version}
      mv $out/lib/perl5/x86_64-linux-thread-multi $out/lib/perl5/site_perl/${perl.version}
    '';
  };

  BioPerl = buildPerlPackage rec {
    name = "BioPerl-1.7.4";
    src = fetchurl {
      url = "mirror://cpan/authors/id/C/CD/CDRAUG/${name}.tar.gz";
      sha256 = "0yvhgifs8g9rwdcq84zw4b005nq2jml6c75zgjscv6d2pd3lj1ss";
    };
    propagatedBuildInputs = with perlPackages; [DBI DataStag Error GD Graph HTTPMessage HTTPMessage IOstringy IOString IPCRun LWP ListMoreUtils SetScalar TestMost TestRequiresInternet URI XMLDOM XMLDOMXPath XMLLibXML XMLLibXML libxml_perl XMLSAX XMLSAXBase XMLSAXWriter XMLTwig XMLWriter YAML DBFile ];
    nativeBuildInputs = with perlPackages; [ TestException TestWarn TestDifferences TestDeep ];
  };

  DataStag = buildPerlPackage rec {
    name = "Data-Stag-0.14";
    src = fetchurl {
      url = "mirror://cpan/authors/id/C/CM/CMUNGALL/${name}.tar.gz";
      sha256 = "0ncf4l39ka23nb01jlm6rzxdb5pqbip01x0m38bnvf1gim825caa";
    };
    propagatedBuildInputs = with perlPackages; [IOString Graph XMLLibXSLT ];
  };

  XMLDOMXPath = buildPerlPackage rec {
    name = "XML-DOM-XPath-0.14";
    src = fetchurl {
      url = "mirror://cpan/authors/id/M/MI/MIROD/${name}.tar.gz";
      sha256 = "1si9m1pqih3ibbd6jnw69fh98dd4krxpx90p65x9j4aja55afwq1";
    };
    propagatedBuildInputs= with perlPackages; [XMLDOM XMLXPathEngine];
    doCheck = false;
  };

  pcapCore = buildPerlPackage rec {
    name = "PCAP-core-${version}";
    version = "4.2.3";
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "PCAP-core";
      rev = version;
      sha256 = "0ia9z0k4jpa02596smivpkw85vi9sv7sbxnysv7jdyh2g1rv1s5c";
    };
    nativeBuildInputs = with perlPackages; [ TestFatal TestWarn ];
    propagatedBuildInputs = with perlPackages; [
      bwa
      samtools
      biobambam2
      cgpBigWig
      ConstFast
      FileWhich
      IPCSystemSimple
      CaptureTiny
      TermUI
      BioDBHTS
      DataUUID
      YAML
      BioPerl
      JSON
      psmisc
    ];
    preConfigure = ''
      sed -i 's|/usr/bin/time|${time}/bin/time|' lib/PCAP/Threaded.pm
      sed -i 's|/bin/bash|${stdenv.shell}|' lib/PCAP/Threaded.pm
    '';
  };

  cgpVcf = buildPerlPackage rec {
    name = "cgpVCF-${version}";
    version = "2.0.4";
    outputs = [ "out" "dev" ];
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "cgpVcf";
      rev = "v${version}";
      sha256 = "0j8spk02v4l0nrqsa42d57djb4nm71daz5hvlm5fwakfq8037yx4";
    };
    nativeBuildInputs = with perlPackages; [ TestFatal TestWarn ];
    propagatedBuildInputs = with perlPackages; [
      ConstFast
      samtools
      vcftools
      DataUUID
      DateTime
      IPCSystemSimple
    ];
  };

  BDebug = buildPerlPackage rec {
    name = "B-Debug-1.26";
    src = fetchurl {
      url = "mirror://cpan/authors/id/R/RU/RURBAN/${name}.tar.gz";
      sha256 = "0fhdaxpkirgnwivd0a7x83dvfwqbngsq2rcfwvfxipgh6i8kyvcd";
    };
  };

  DevelCover = buildPerlPackage rec {
    name = "Devel-Cover-1.31";
    src = fetchurl {
      url = "mirror://cpan/authors/id/P/PJ/PJCJ/${name}.tar.gz";
      sha256 = "0wsd7fkpy5qwnz68hz624frnbqii5igg6l2q7n1mx1lh32vn0kcr";
    };
    nativeBuildInputs = with perlPackages; [ BDebug ];
    doCheck = false;
  };


  alleleCount = stdenv.mkDerivation rec {
    name = "alleleCount-${version}";
    version = "4.0.1";
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "alleleCount";
      rev = "v${version}";
      sha256 = "0nkwnjqglgshzhlmz1r0khdjai9mfz4ih8bzrzg0g18d1725k6gp";
    };
    buildInputs = [ htslib perl zlib bzip2 lzma ];
    preConfigure = "cd c";
    installPhase = ''
      mkdir -p $out/bin
      cp bin/* $out/bin
    '';
  };

  alleleCountPerl = buildPerlPackage rec {
    name = "alleleCount.pl-${version}";
    version = "4.0.1";
    src = fetchFromGitHub {
    owner = "cancerit";
    repo = "alleleCount";
    rev = "v${version}";
    sha256 = "0nkwnjqglgshzhlmz1r0khdjai9mfz4ih8bzrzg0g18d1725k6gp";
    };
    preConfigure = ''
      cd perl
    '';
    nativeBuildInputs = with perlPackages; [
      TestFatal
    ];
    propagatedBuildInputs = with perlPackages; [
      ConstFast
      DevelCover
      FileSlurp
      FileWhich
      PodCoverage
      TryTiny
      alleleCount
      IPCSystemSimple
    ];
  };

  ascatR = rWrapper.overrideAttrs (attrs: {
    packages = with rPackages; [ RColorBrewer ];
  });

  ascatNGS = buildPerlPackage rec {
    name = "ascatNGS-${version}";
    version = "4.1.2";
    src = fetchFromGitHub {
      owner = "cancerit";
      repo = "ascatNgs";
      rev = "v${version}";
      sha256 = "0qc1wp37k2c9vya992cfiphp2jpp8hvfvqz5ydi2d5b71prfaw82";
    };
    preConfigure = ''
      cd perl
      cp ${ascat} share/ascat/ascat.R
    '';
    nativeBuildInputs = with perlPackages; [ TestFatal TestWarn ];
    propagatedBuildInputs = with perlPackages; [
      ascatR
      FileShareDirInstall
      pcapCore
      FileShareDir
      cgpVcf
      alleleCountPerl
    ];
  };

in ascatNGS
