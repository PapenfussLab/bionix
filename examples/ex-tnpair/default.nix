{ bionix ? import <bionix> { } }:

with bionix;
with pkgs;
with lib;

let
  pair = {
    normal = {
      type = "reference";
      inputs = {
        input1 = {
          url =
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2752449/ERR2752449_1.fastq.gz";
          sha256 =
            "52f8b1b1a58b60c66ce566371dfe7a1301a787e8521a4ee41019bbf4f4d18dfe";
        };
        input2 = {
          url =
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2752449/ERR2752449_2.fastq.gz";
          sha256 =
            "9d1e2ea772bbdf5ff3ee6a44d2d4244155b7d195a37745a2028628e2543cd8f0";
        };
      };
    };

    tumour = {
      type = "melanoma";
      inputs = {
        input1 = {
          url =
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2752450/ERR2752450_1.fastq.gz";
          sha256 =
            "2b3c98c36c2b2b6bc4682401a592a900f8eb2a143f93494ee448d6b075c12ec7";
        };
        input2 = {
          url =
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2752450/ERR2752450_2.fastq.gz";
          sha256 =
            "0569beded708ef520dadca45ab8a70bd890caf441a0ad3749397f315dc1d2e8c";
        };
      };
    };
  };

  fetch = s: mapAttrs (_: fetchFastQGZ) s.inputs;

in
import ./tnpair.nix {
  inherit pair fetch bionix;
}
