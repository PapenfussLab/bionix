{bionix, nixpkgs}:

with nixpkgs;
with bionix;

{
  uncompress = f: types.matchFiletype "uncompress" {
    fa = _: f;
    fq = _: f;
    bam = _: f;
    sam = _: f;
    cram = _: f;
    vcf = _: f;
    bed = _: f;
    gz = _: types.tagFiletype (types.gunzip f.filetype) (stdenv.mkDerivation {
      name = "gunzip";
      buildCommand = "gunzip < ${f} > $out";
    });
    bz2 = _: types.tagFiletype (types.bunzip2 f.filetype) (stdenv.mkDerivation {
      name = "bunzip2";
      buildCommand = "bunzip2 < ${f} > $out";
    });
  } f.filetype;

  gzip = f:
    let
      gz = (stdenv.mkDerivation {
        name = "gzip";
        buildCommand = "gzip < ${f} > $out";
        passthru = { filetype = types.filetype.gz f.filetype; };
      });
    in types.matchFiletype "compressed" {
      fa = _: gz;
      fq = _: gz;
      bam = _: gz;
      sam = _: gz;
      cram = _: gz;
      vcf = _: gz;
      bed = _: gz;
    } f;

  bzip2 = f:
    let
      bz2 = (stdenv.mkDerivation {
        name = "bzip2";
        buildCommand = "bzip2 < ${f} > $out";
        passthru = { filetype = types.filetype.bz2 f.filetype; };
      });
    in types.matchFiletype "compressed" {
      fa = _: gz;
      fq = _: gz;
      bam = _: gz;
      sam = _: gz;
      cram = _: gz;
      vcf = _: gz;
      bed = _: gz;
    } f;
}
