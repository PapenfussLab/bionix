# This example uses the pipelines specified in the call.nix file on the
# synthetic data in this directory.
with import <bionix> {};

let

  # List of input samples
  inputs = [
    # Sample 1
    {
      input1 = fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/example/sample1-1.fq";
        sha256 = "1m3vc248mbr4v56459q4xsklznssgqb35lwhwk1i9qvxqm7c72nq";
      };

      input2 =  fetchFastQ {
      url = "https://github.com/PapenfussLab/bionix/raw/master/example/sample1-2.fq";
        sha256 = "13fqvdb5r2sidi2i0s3ifg8gyxp8kibpxc3cbiw07d5zcn3g479x";
      };
    }

    # Sample 2
    {
      input1 = fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/example/sample2-1.fq";
        sha256 = "0c6ik4smsw2kb8xbbci80d00c24klk6mqd2a71y4v9hmpnc42fmr";
      };

      input2 =  fetchFastQ {
        url = "https://github.com/PapenfussLab/bionix/raw/master/example/sample2-2.fq";
        sha256 = "0zvrkm2m69lwfi89hzxsbnkp4i6cszh7c624l1qwmd4s6gh7vfcx";
      };
    }
  ];

  # The reference for the synthetic data
  ref = fetchFastA {
    url = "https://github.com/PapenfussLab/bionix/raw/master/example/ref.fa";
    sha256 = "06gphhh40h3mvwvs2m51qc3rpih8mcs8frhd48l94d5bzwfhb2hc";
  };

in import ./call.nix {inherit inputs ref;}
