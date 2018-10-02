{ bionix, nixpkgs }:

with nixpkgs;

rec {
  grch38 = grch38-p12;
  grch38-p12 = {
    seq = stdenvNoCC.mkDerivation rec {
        name = "seq-grch38.${version}";
        version = "p12";
        src = fetchurl {
          url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz";
          sha256 = "0ji2ggpmgnbpwbhq8mirj6h3lyy02nl2rnz7n892iq5cqpsblh4z";
        };
        buildCommand = "gunzip < $src > $out";
      };
    blacklist = stdenvNoCC.mkDerivation {
      name = "blacklist-grch38";
      src = fetchurl {
        url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz";
        sha256 = "1lpnqq1mjidbdxj5i6x26rxa8x1rs8q3hlf0z1z49j3jsnkgffky";
      };
      buildCommand = "gunzip < $src > $out";
    };
  };

  mm10 = mm10-p4;
  mm10-p4 = {
      seq = stdenvNoCC.mkDerivation rec {
          name = "seq-mm10.${version}";
          version = "p4";
          src = fetchurl {
              url = "ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/mm10Patch4/mm10Patch4.fa.gz";
              md5 = "0dfc06235d8e96b8e364932c8fc02456";
              };
          buildCommand = "gunzip < $src > $out";
      };
  };
}
