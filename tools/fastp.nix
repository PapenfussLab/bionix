{ bionix }:

with bionix;

rec {
    app = pkgs.callPackage ./fastp-app.nix {};
    check = callBionixE ./fastp-check.nix;
}

