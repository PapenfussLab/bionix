{ bionix, nixpkgs }:

with nixpkgs;

rec {
  grch38 = grch38-p12;
  grch38-p12 = {
    seq = stdenv.mkDerivation rec {
        name = "seq-grch38.${version}";
        version = "p12";
        src = fetchurl {
          url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz";
          sha256 = "0ji2ggpmgnbpwbhq8mirj6h3lyy02nl2rnz7n892iq5cqpsblh4z";
        };
        buildCommand = "gunzip < $src > $out";
      };
  };
}
