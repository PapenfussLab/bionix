{ bionix }:

with bionix;
with types;

let
  gzip = callBionixE ./compression-gzip.nix;
  gunzip = callBionixE ./compression-gunzip.nix;
  bunzip2 = callBionixE ./compression-bunzip2.nix;

in
{
  uncompress = attrs: f: matchFiletype "uncompress"
    {
      fa = _: f;
      fq = _: f;
      bam = _: f;
      sam = _: f;
      cram = _: f;
      vcf = _: f;
      bed = _: f;
      gz = _: gunzip attrs f;
      bgz = _: gunzip attrs f;
      bz2 = _: bunzip2 attrs f;
    }
    f;

  gzip = attrs: f:
    let gz = gzip attrs f;
    in
    types.matchFiletype "compressed"
      {
        fa = _: gz;
        fq = _: gz;
        bam = _: gz;
        sam = _: gz;
        cram = _: gz;
        vcf = _: gz;
        bed = _: gz;
        gz = _: f;
        bgz = _: f;
        bz2 = _: with compression; gzip { } (uncompress { } f);
      }
      f;

  bgzip = attrs: f:
    let gz = gzip (attrs // { block = true; }) f;
    in
    types.matchFiletype "compressed"
      {
        fa = _: gz;
        fq = _: gz;
        bam = _: gz;
        sam = _: gz;
        cram = _: gz;
        vcf = _: gz;
        bed = _: gz;
        bgz = _: f;
        gz = _: with compression; bgzip { } (uncompress { } f);
        bz2 = _: with compression; bgzip { } (uncompress { } f);
      }
      f;


  bzip2 = attrs: f:
    let bz2 = bzip2 attrs f;
    in
    types.matchFiletype "compressed"
      {
        fa = _: gz;
        fq = _: gz;
        bam = _: gz;
        sam = _: gz;
        cram = _: gz;
        vcf = _: gz;
        bed = _: gz;
        bz2 = _: f;
        bgz = _: with compression; bzip2 { } (uncompress { } f);
        gz = _: with compression; bzip2 { } (uncompress { } f);
      }
      f;
}
