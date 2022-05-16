{bionix}:
with bionix;

{
app = callPackage ./quip-app.nix {};
quip = callBionix ./quip-quip.nix;
unquip = callBionix ./quip-unqip.nix;
}
