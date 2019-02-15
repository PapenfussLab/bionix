{ bionix }:

with bionix;
with types;

rec {
  grch38 = grch38-p12;
  grch38-p12 = rec {
    seq = stage rec {
      name = "seq-grch38.${version}";
      version = "p12";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz";
        sha256 = "0ji2ggpmgnbpwbhq8mirj6h3lyy02nl2rnz7n892iq5cqpsblh4z";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.fa {};
    };
    blacklist = stage {
      name = "blacklist-grch38";
      src = pkgs.fetchurl {
        url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz";
        sha256 = "1lpnqq1mjidbdxj5i6x26rxa8x1rs8q3hlf0z1z49j3jsnkgffky";
      };
      buildCommand = "gunzip < $src > $out";
      passthru.filetype = filetype.bed { ref = seq; };
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
      cosmic = {coding, noncoding}: stage rec {
        name = "cosmic-grch38";
        buildInputs = with pkgs; [ gawk ];
        buildCommand = ''
        gunzip < ${coding} | grep '^#' > $out
        cat ${coding} ${noncoding} | gunzip | grep -v '^#' | sed 's/^/chr/' | sort -t$'\t' -k1,1 -k2,2n >> $out
        '';
        passthru.filetype = filetype.vcf { ref = seq; };
      };
      ensembl = {
        cdna = stage rec {
          name = "ensembl-grch38-cdna-${version}";
          version = "94";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";
            sha256 = "1fc5d6p2wlwsm49wnmxmm3byjx5jvr6z9fpzrq7v7fpb086adl0h";
          };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
        };
        ncrna = stage rec {
          name = "ensembl-grch38-ncrna-${version}";
          version = "94";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz";
            sha256 = "1cpasykwriila52nqgvw6d3mjyh6d9qi613hvhn4h1dxkqzgnjff";
          };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
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
          db = pkgs.requireFile {
            name = "dbNSFP.txt.gz";
            message = "download the dbNSFP database manually from https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlbTZodjlGUDZnTGc and add to nix store";
            sha256 = "0gahnwkc7v2q6p6ixkhvsgqvvm6xf0c3bdh4nf0alih83h3wffd0";
          };
          index = pkgs.requireFile {
            name = "dbNSFP.txt.gz.tbi";
            message = "download the dbNSFP index manually from https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlNVBJdFA5cFZRYkE and add to nix store";
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
      ensembl = {
        cdna = stage rec {
          name = "ensembl-grch38-cdna-${version}";
          version = "94";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz";
            sha256 = "0khp9l6s35lav2xqp7vkk6ybnz4wjihn7lapjf2lbpnbzjb4hp6d";
          };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
        };
        ncrna = stage rec {
          name = "ensembl-grch38-ncrna-${version}";
          version = "94";
          src = pkgs.fetchurl {
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
