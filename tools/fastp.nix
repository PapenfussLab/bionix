{ bionix }:

with bionix;

rec {
  /* Check and filter fastqs
    Type: { ... } -> { input1, input2 ? null } -> { out :: html, fastq1 :: fastq, fastq2 :: fastq, json :: JSON } 
  */
  check = callBionixE ./fastp-check.nix;
}

