{bionix}:

with bionix;
with lib;

let
  nix-adt-src = pkgs.fetchFromGitHub {
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

  defError = errF: y: x: listToAttrs (map (n: nameValuePair n (errF n)) (filter (x: builtins.substring 0 1 x != "_") (attrNames x))) // y;

in
rec {
  nix-adt = import "${nix-adt-src}";
  matchFiletype = sym: y: x: if x ? filetype then matchFiletype' sym y x.filetype else abort "unknown filetype for ${sym}";
  matchFiletype' = sym: y: x: match x (defError (idft sym) y filetype);
  filetype = make-type "filetype" {
    fa = {};
    fq = {};
    bam = {ref = any; sorting = sort;};
    sam = {ref = any; sorting = sort;};
    cram = {ref = any; sorting = sort;};
    vcf = {ref = any;};
    bed = {ref = any;};
    gz = filetype;
    bz2 = filetype;
  };

  toCram = matchFiletype "bam2cram" { bam = filetype.cram; sam = filetype.cram; cram = filetype.cram; };
  toBam = matchFiletype "bam2cram" { bam = filetype.bam; sam = filetype.bam; cram = filetype.bam; };
  toSam = matchFiletype "bam2cram" { bam = filetype.sam; sam = filetype.sam; cram = filetype.sam; };

  matchSorting = sym: y: x: match x.sorting (defError (idst sym) y sort);
  matchFileSorting = sym: y: let f = matchSorting sym y; in matchFiletype sym { bam = f; sam = f; cram = f; };
  sort = make-type "sort" {
    none = {};
    coord = {};
    name = {};
  };
  coordSort = f: matchFiletype "coordSort" { bam = x: filetype.bam (x // {sorting = sort.coord {};}); } {filetype = f;};
  nameSort = f: matchFiletype "nameSort" { bam = x: filetype.bam (x // {sorting = sort.name {};}); } {filetype = f;};

  gunzip = matchFiletype "gunzip" { gz = x: x; };
  bunzip2 = matchFiletype "bunzip2" { bz2 = x: x; };

  tag = attrs: x: if x ? type && x.type == "derivation" then x // attrs else tagPassthru attrs x;
  tagPassthru = attrs: x: if x ? passthru then x // { passthru = x.passthru // attrs; } else x // { passthru = attrs; };
  tagFiletype = ft: tag { filetype = ft; };
}
