{ bionix }:

with bionix;

rec {
    app = pkgs.callPackage ./fastp-app.nix {};

    /* Check and filter fastqs
    Type: { ... } -> { input1, input2 ? null } -> { out :: html, fastq1 :: fastq, fastq2 :: fastq, json :: JSON } 
    */
    check = callBionixE ./fastp-check.nix;
}

