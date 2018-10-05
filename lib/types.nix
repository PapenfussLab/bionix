{bionix, nixpkgs}:

with nixpkgs;

let
  nix-adt-src = fetchFromGitHub {
    owner = "shlevy";
    repo = "nix-adt";
    rev = "dd04b5d08eed65ecd73aafde56a78078e09f1c74";
    sha256 = "0vhk1y7gilgn2pgvj4813bh2ljpw4pvrph8k8b0fsg56dbm8mqxa";
  };
  nix-adt = import "${nix-adt-src}";
  inherit (nix-adt.checked) make-type match any std none;
  inherit (std) option;

  idft = sym: ft: _: abort "unhandled filetype (${ft}) for ${sym}";
  idst = sym: st: _: abort "unhandled sorting (${st}) for ${sym}";

in
rec {
  option-sort = option sorting;

  matchFiletype = sym: y: x: if x ? filetype then match x.filetype ({
    fa = idft sym "fasta";
    fq = idft sym "fastq";
    bam = idft sym "bam";
    sam = idft sym "sam";
    cram = idft sym "cram";
    vcf = idft sym "vcf";
    bed = idft sym "bed";
    gz = idft sym "gz";
    bz2 = idft sym "bz2";
  } // y) else abort "unknown filetype for ${sym}";
  filetype = make-type "filetype" {
    fa = {};
    fq = {};
    bam = {ref = any; sorting = option-sort;};
    sam = {ref = any; sorting = option-sort;};
    cram = {ref = any; sorting = option-sort;};
    vcf = {ref = any;};
    bed = {ref = any;};
    gz = filetype;
    bz2 = filetype;
  };

  toCram = matchFiletype "bam2cram" { bam = filetype.cram; sam = filetype.cram; cram = filetype.cram; };
  toBam = matchFiletype "bam2cram" { bam = filetype.bam; sam = filetype.bam; cram = filetype.bam; };
  toSam = matchFiletype "bam2cram" { bam = filetype.sam; sam = filetype.sam; cram = filetype.sam; };

  matchSorting = sym: y: let f = x: match x.sorting { some = z: match z ( { coord = idst sym "coord"; name = idst sym "name"; } // y); none = abort "unknown sort for ${sym}"; }; in matchFiletype sym { bam = f; sam = f; cram = f; };
  sorting = make-type "sorting" {
    coord = {};
    name = {};
  };
  coordSort = f: matchFiletype "coordSort" { bam = x: filetype.bam (x // {sorting = option-sort.some (sorting.coord {});}); } {filetype = f;};
  nameSort = f: matchFiletype "nameSort" { bam = x: filetype.bam (x // {sorting = option-sort.some (sorting.name {});}); } {filetype = f;};

  gunzip = matchFiletype "gunzip" { gz = x: x; };
  bunzip2 = matchFiletype "bunzip2" { bz2 = x: x; };

  tag = attrs: x: if x ? type && x.type == "derivation" then x // attrs else tagPassthru attrs x;
  tagPassthru = attrs: x: if x ? passthru then x // { passthru = x.passthru // attrs; } else x // { passthru = attrs; };
  tagFiletype = ft: tag { filetype = ft; };
}
