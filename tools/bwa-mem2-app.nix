{ bwa, fetchFromGitHub }:
  bwa.overrideAttrs (_:
    {
      name = "bwa-mem2";
      src = fetchFromGitHub {
        owner = "bwa-mem2";
        repo = "bwa-mem2";
        rev = "f882015f7f46845f72c10094a621264f78d206d8";
        sha256 = "1wfx8j9mwbb29jw4zxp28lajmj774jy1l5x229nf065cqr2dgyqg";
        };
      installPhase = ''
        mkdir -p $out/bin
        cp bwa-mem2 $out/bin
        '';
      })
