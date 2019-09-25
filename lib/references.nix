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
        name = "GRCh37.87";
        src = pkgs.fetchurl {
          url = "mirror://sourceforge/project/snpeff/databases/v4_3/snpEff_v4_3_${name}.zip";
          sha256 = "0ybbj4470ilc4csmgfjqd6hqq4krwjws97ywjnqrhbi4dcq3h3bg";
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
          message = "download the dbNSFP database manually from https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlSTY5dDJjcHVRZ3M and add to nix store";
          sha256 = "0gfzbid3pc10zds7ya50w4qfynsxgpyh7dx35vhm5f3h64mw75pm";
        };
        index = pkgs.requireFile {
          name = "dbNSFP.txt.gz.tbi";
          message = "download the dbNSFP index manually from https://drive.google.com/uc?export=download&id=0B7Ms5xMSFMYlOTV5RllpRjNHU2s and add to nix store";
          sha256 = "0bwwigbnz32mmc8bczidjf68vv8x8i28zwkl2kgcbpj542zk5q86";
        };
      };
    };
  };
  grch38 = grch38-p12;
  grch38-p12 = rec {
    seq = stage rec {
      name = "seq-grch38.${version}";
      version = "32";
      src = pkgs.fetchurl {
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}/GRCh38.p13.genome.fa.gz";
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
      ensembl = let version = "94"; in {
        cdna = stage {
          name = "ensembl-grch38-cdna-${version}";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";
            sha256 = "1fc5d6p2wlwsm49wnmxmm3byjx5jvr6z9fpzrq7v7fpb086adl0h";
          };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
        };
        ncrna = stage {
          name = "ensembl-grch38-ncrna-${version}";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz";
            sha256 = "1cpasykwriila52nqgvw6d3mjyh6d9qi613hvhn4h1dxkqzgnjff";
          };
          buildCommand = "gunzip < $src > $out";
          passthru.filetype = filetype.fa {};
        };
        gtf = stage {
          name = "ensembl-grch38-gtf-${version}";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz";
            sha256 = "124swjp6bhhl0wjhfq16675s3fqyny293mhrjwvzkiybpgfwzpi7";
          };
          buildCommand = "gunzip < $src > $out";
        };
        gff3 = stage {
          name = "ensembl-grch38-gff-${version}";
          src = pkgs.fetchurl {
            url = "ftp://ftp.ensembl.org/pub/release-${version}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${version}.gff3.gz";
            sha256 = "0r2npkpq7z0dd74fl8q3kzv0860chfx3v319z6libl5qy9rspdmi";
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
