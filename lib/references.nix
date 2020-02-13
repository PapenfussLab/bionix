{ bionix }:

with bionix;
with types;

rec {
  grch37 = rec {
    seq = stage rec {
      name = "seq-grch37.${version}";
      version = "19";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/GRCh37.p13.genome.fa.gz";
        sha256 = "1midmq3kaci3nizg37nqzv9syfasimgxxvg8az8gz0fi1ksljw53";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    dbsnp = stage {
      name = "dbsnp-b151_GRCh37p13";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz";
        sha256 = "0jf71h5wy82xhgsjkvp05mj2grjrjlnswmr0wz4lb87g3ip3c2mm";
      };
      buildInputs = with pkgs; [ gawk ];
      buildCommand = ''
      gunzip < $src | awk '/^[^#]/{print "chr" $0;next}{print}' > $out
        '';
        passthru.filetype = filetype.vcf { ref = seq; };
      };
    ensembl = let version = "74"; in {
      cdna = stage {
        name = "ensembl-grch37-cdna-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.${version}.cdna.all.fa.gz";
          sha256 = "1m62hiw17zcxg3sza0aq53885wb8g202j8lc1ilhmkg2izzbyihj";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      gtf = stage {
        name = "ensembl-grch37-gtf-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh37.${version}.gtf.gz";
          sha256 = "1m62hiw17zcxg3sza0aq53885wb8g202j8lc1ilhmkg2izzbyihj";
        };
        buildCommand = "gunzip < $src > $out";
      };
    };
    snpeff = {
      db = pkgs.stdenv.mkDerivation rec {
        name = "GRCh37.75";
        src = pkgs.fetchurl {
          url = "mirror://sourceforge/project/snpeff/databases/v4_3/snpEff_v4_3_${name}.zip";
          sha256 = "19c8wwx91vq47z7j7f455vsv8jw067x5rd7449d1z0nln82zpmhm";
        };
        buildInputs = with pkgs; [ unzip ];
        buildCommand = ''
          unzip ${src}
          mv data/${name} $out
        '';
      };
      dbnsfp = {
        db = fetchgdrive {
          url = "https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlSTY5dDJjcHVRZ3M";
          sha256 = "0gfzbid3pc10zds7ya50w4qfynsxgpyh7dx35vhm5f3h64mw75pm";
        };
        index = fetchgdrive {
          url = "https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlOTV5RllpRjNHU2s";
          sha256 = "0bwwigbnz32mmc8bczidjf68vv8x8i28zwkl2kgcbpj542zk5q86";
        };
      };
    };
  };
  grch38 = rec {
    seq = stage rec {
      name = "seq-grch38.${version}";
      version = "31";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/GRCh38.p12.genome.fa.gz";
        sha256 = "1cmpg634likx3gbwn578491jjb8k7xm3cmqaqs3r9hyd195vbw98";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    dbsnp = stage {
      name = "dbsnp-b151_GRCh38p7";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz";
        sha256 = "0r6m2yrcfw8bbdca515axjls30ssjas6x3qwi5qz07l3prjwmdd4";
      };
      buildInputs = with pkgs; [ gawk ];
      buildCommand = ''
      gunzip < $src | awk '/^[^#]/{print "chr" $0;next}{print}' > $out
        '';
        passthru.filetype = filetype.vcf { ref = seq; };
      };
    ensembl = let version = "97"; in {
      cdna = stage {
        name = "ensembl-grch38-cdna-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";
          sha256 = "1m6hvnvlrsi6bzcmq0lnv0igy3in1a7jp723yc74g4g6zjp3cy8c";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      ncrna = stage {
        name = "ensembl-grch38-ncrna-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz";
          sha256 = "1r0dmybn31wf6xc90z2c08ngivkv39hqa8wqg3vik6s4spwpdhj0";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      gtf = stage {
        name = "ensembl-grch38-gtf-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz";
          sha256 = "068ab5jf87il301jcr0576c4q0smv6kxpv94gnrm3qzl6kvmaawd";
        };
        buildCommand = "gunzip < $src > $out";
      };
      gff3 = stage {
        name = "ensembl-grch38-gff-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${version}.gff3.gz";
          sha256 = "1xvhsn938mw0032qgc9dvw3k2xrhpx77b8ms03fkrs2s67f7zli7";
        };
        buildCommand = "gunzip < $src > $out";
      };
    };
    snpeff = {
      db = pkgs.stdenv.mkDerivation rec {
        name = "GRCh38.86";
        src = pkgs.fetchurl {
          url = "mirror://sourceforge/project/snpeff/databases/v4_3/snpEff_v4_3_${name}.zip";
          sha256 = "1rf8q7l732ayjq2lpny4s75zpij05j00151374nqblk4wri2mz0i";
        };
        buildInputs = with pkgs; [ unzip ];
        buildCommand = ''
          unzip ${src}
          mv data/${name} $out
        '';
      };
      dbnsfp = {
        db = fetchgdrive {
          url = "https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlbTZodjlGUDZnTGc";
          sha256 = "0gahnwkc7v2q6p6ixkhvsgqvvm6xf0c3bdh4nf0alih83h3wffd0";
        };
        index = fetchgdrive {
          url = "https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlNVBJdFA5cFZRYkE";
          sha256 = "18blkly6gvg7r0sx968xlb1zl2kqg5j1kpbrm2r7ajlxlfyvrx3w";
        };
      };
    };
    UCSCgenes = stage {
      name = "UCSCgenes";
      src = pkgs.fetchurl {
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz";
        sha256 = "1f53myn4vpvswzssx2xsiq9si8w58gpcm0f32srq220w36hq9md4";
      };
      buildCommand = "gunzip < $src > $out";
    };
  };

  grcm38 = grcm38-p6;
  grcm38-p6 = {
    seq = stage rec {
      name = "seq-grcm38.${version}";
      version = "p6";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.${version}.genome.fa.gz";
        sha256 = "0ryiqab5bldpzawylsk2qpjxr2j701q03ww9jqyxhkimqpn9g3mr";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    ensembl = let
      version = "94";
    in {
      cdna = stage {
        name = "ensembl-grcm38-cdna-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz";
          sha256 = "0khp9l6s35lav2xqp7vkk6ybnz4wjihn7lapjf2lbpnbzjb4hp6d";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      ncrna = stage {
        name = "ensembl-grcm38-ncrna-${version}";
        version = "94";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz";
          sha256 = "0d997gm8p2b89rm5d46m2x4vz9lijxarfr2lzylnbi8gyqrbagdd";
        };
        buildCommand = "gunzip < $src > $out";
        passthru.filetype = filetype.fa {};
      };
      gtf = stage {
        name = "ensembl-grcm38-gtf-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/gtf/mus_musculus/Mus_musculus.GRCm38.${version}.gtf.gz";
          sha256 = "0i61jq5i5bcini5nxqxxp3rnz2xzgychvzdn0k451f5rv053lp3v";
        };
        buildCommand = "gunzip < $src > $out";
      };
      gff3 = stage {
        name = "ensembl-grcm38-gff-${version}";
        src = pkgs.fetchurl {
          url = "ftp://ftp.ensembl.org/pub/release-${version}/gff3/mus_musculus/Mus_musculus.GRCm38.${version}.gff3.gz";
          sha256 = "15fmdpx6g96fygwhs10jwrb2q5p9y64bc3d4clg856k57qzzgprg";
        };
        buildCommand = "gunzip < $src > $out";
      };
    };
  };

  mm10 = mm10-p4;
  mm10-p4 = {
    seq = stage rec {
      name = "seq-mm10.${version}";
      version = "p4";
      src = pkgs.fetchurl {
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/mm10Patch4/mm10Patch4.fa.gz";
        sha256 = "1660d6d05f3aa266c6053cfd1efef1747d9e854836917241d6f47cff7a55340c";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
  };
}
