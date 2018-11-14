{ bionix, nixpkgs }:

with nixpkgs;
with bionix.types;

rec {
  grch38 = grch38-p12;
  grch38-p12 = rec {
    seq = stdenvNoCC.mkDerivation rec {
      name = "seq-grch38.${version}";
      version = "p12";
      src = fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz";
        sha256 = "0ji2ggpmgnbpwbhq8mirj6h3lyy02nl2rnz7n892iq5cqpsblh4z";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    blacklist = stdenvNoCC.mkDerivation {
      name = "blacklist-grch38";
      src = fetchurl {
        url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz";
        sha256 = "1lpnqq1mjidbdxj5i6x26rxa8x1rs8q3hlf0z1z49j3jsnkgffky";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.bed { ref = seq; };
    };
    dbsnp = stdenvNoCC.mkDerivation {
      name = "dbsnp-b151_GRCh38p7";
      src = fetchurl {
        url = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz";
        sha256 = "0r6m2yrcfw8bbdca515axjls30ssjas6x3qwi5qz07l3prjwmdd4";
      };
      buildInputs = [ gawk ];
      buildCommand = ''
        gunzip < $src | awk '/^[^#]/{print "chr" $0;next}{print}' > $out
      '';
      passthru.filetype = filetype.vcf { ref = seq; };
    };
    cosmic = {coding, noncoding}: stdenvNoCC.mkDerivation rec {
      name = "cosmic-grch38";
      buildInputs = [ gawk ];
      buildCommand = ''
        gunzip < ${coding} | grep '^#' > $out
        cat ${coding} ${noncoding} | gunzip | grep -v '^#' | sed 's/^/chr/' | sort -t$'\t' -k1,1 -k2,2n >> $out
      '';
      passthru.filetype = filetype.vcf { ref = seq; };
    };
    ensembl = {
      cdna = stdenvNoCC.mkDerivation rec {
        name = "ensembl-grch38-cdna-${version}";
        version = "94";
        src = fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";
          sha256 = "1fc5d6p2wlwsm49wnmxmm3byjx5jvr6z9fpzrq7v7fpb086adl0h";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      ncrna = stdenvNoCC.mkDerivation rec {
        name = "ensembl-grch38-ncrna-${version}";
        version = "94";
        src = fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz";
          sha256 = "1cpasykwriila52nqgvw6d3mjyh6d9qi613hvhn4h1dxkqzgnjff";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
    };
  };

  grcm38 = grcm38-p6;
  grcm38-p6 = {
    seq = stdenvNoCC.mkDerivation rec {
      name = "seq-grcm38.${version}";
      version = "p6";
      src = fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.${version}.genome.fa.gz";
        sha256 = "0ryiqab5bldpzawylsk2qpjxr2j701q03ww9jqyxhkimqpn9g3mr";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    ensembl = {
      cdna = stdenvNoCC.mkDerivation rec {
        name = "ensembl-grch38-cdna-${version}";
        version = "94";
        src = fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz";
          sha256 = "0khp9l6s35lav2xqp7vkk6ybnz4wjihn7lapjf2lbpnbzjb4hp6d";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      ncrna = stdenvNoCC.mkDerivation rec {
        name = "ensembl-grch38-ncrna-${version}";
        version = "94";
        src = fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz";
          sha256 = "0d997gm8p2b89rm5d46m2x4vz9lijxarfr2lzylnbi8gyqrbagdd";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
    };
  };

  mm10 = mm10-p4;
  mm10-p4 = {
      seq = stdenvNoCC.mkDerivation rec {
          name = "seq-mm10.${version}";
          version = "p4";
          src = fetchurl {
              url = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/mm10Patch4/mm10Patch4.fa.gz";
              sha256 = "1660d6d05f3aa266c6053cfd1efef1747d9e854836917241d6f47cff7a55340c";
              };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
      };
  };
}
