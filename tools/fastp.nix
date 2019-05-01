{ bionix }:

with bionix;

rec {
    app = pkgs.callPackage ./fastp-app.nix {};
    run = callBionixE ./fastp-run.nix;
}

