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
}
