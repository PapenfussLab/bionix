{ bionix, normalName, genomeVersion ? "hg38" }:
vcf:

with bionix;
with pkgs;

let
  refMap = with rPackages; {
    hg38 = BSgenome_Hsapiens_UCSC_hg38;
    hg19 = BSgenome_Hsapiens_UCSC_hg19;
  };

  renv = rWrapper.override {
    packages = with rPackages; [
      (refMap."${genomeVersion}" or (abort "unsupported reference"))
      Biobase
      BiocGenerics
      Biostrings
      DelayedArray
      IRanges
      S4Vectors
      StructuralVariantAnnotation
      VariantAnnotation
      XVector
      argparser
      ggplot2
      matrixStats
      readr
      rtracklayer
      stringdist
      stringr
      testthat
      tidyverse
    ];
  };

  script = stdenvNoCC.mkDerivation {
    name = "gridss-somatic-filtering-scripts";
    src = fetchFromGitHub {
      owner = "PapenfussLab";
      repo = "gridss";
      rev = "a6b230a78179869c0210dc878490811be813d2fb";
      sha256 = "sha256-b+b6BEvjZbKYR+pJ/Z8yQpcZpairpSUoMHtRRTuRUls=";

    };

    doBuild = false;
    nativeBuildInputs = [ makeWrapper ];

    installPhase = ''
      mkdir -p $out/libexec/gridss
      cp scripts/{gridss.config.R,libgridss.R,gridss_somatic_filter} $out/libexec/gridss
      mkdir -p $out/bin
      makeWrapper ${renv}/bin/Rscript $out/bin/gridss_somatic_filter \
        --add-flags $out/libexec/gridss/gridss_somatic_filter \
        --add-flags "--ref BSgenome.Hsapiens.UCSC.${genomeVersion}" \
        --add-flags "--scriptdir $out/libexec/gridss"
    '';
  };

  findNormal = writeText "find-normal.awk" ''
    /^#C/{
      for(i = 10; i <= NF && $i != "${normalName}"; i++);
      if(i > NF){
        print "findNormal: could not match name" > /dev/stderr
        exit(1)
      }
      printf("%s", i - 9)
      exit(0)
    }
  '';

in
stage {
  name = "gridss-somatic";
  buildCommand = ''
    NORM=$(awk -f ${findNormal} ${vcf})
    ln -s ${vcf} in.vcf
    ${script}/bin/gridss_somatic_filter \
      --output ./out \
      --input in.vcf \
      --normalordinal $NORM
    gunzip < out.bgz > $out
  '';
}
