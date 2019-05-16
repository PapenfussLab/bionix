# This example uses the pipelines specified in the call.nix file on the
# synthetic data in this directory.
{bionix ? import ./.. {}}:

with bionix;

let

  # List of input samples
  inputs = [
    # Sample 1
    {
      input1 = fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/examples/sample1-1.fq";
        sha256 = "0kh29i6fg14dn0fb1xj6pkpk6d83y7zg7aphkbvjrhm82braqkm8";
      };

      input2 =  fetchFastQ {
      url = "https://github.com/PapenfussLab/bionix/raw/master/examples/sample1-2.fq";
        sha256 = "0czk85km6a91y0fn4b7f9q7ps19b5jf7jzwbly4sgznps7ir2kdk";
      };
    }

    # Sample 2
    {
      input1 = fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/examples/sample2-1.fq";
        sha256 = "08gixavfklqvk1m2ic6v56z82vl00qnpsd9xb64z6zl03nz98mcy";
      };

      input2 =  fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/examples/sample2-2.fq";
        sha256 = "1xxwm2vq52axpdhm14rh5mg5nzzpxaqnvhzrqhajm27fqksgzjjw";
      };
    }
  ];

  # The reference for the synthetic data
  ref = fetchFastA {
    url = "https://github.com/PapenfussLab/bionix/raw/master/examples/ref.fa";
    sha256 = "0sy9hq8n55knfkiblam50dzaiwhrx6pv8b8l1njdn6kfj4wflz2p";
  };

in import ./call.nix {inherit inputs ref bionix;}
