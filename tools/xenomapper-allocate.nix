{ bionix
, flags ? null
}:

{primary, secondary}:

with bionix;
with lib;
with types;

let

  isSortedBam = matchFiletype "xenomapper-allocate" {bam = matchSorting "xenomapper-allocate" { coord = _: false; name = _: true; none = _: false; }; };
  outs = [ "primary_specific" "primary_multi" "secondary_specific" "secondary_multi" "unassigned" "unresolved"];

in

assert isSortedBam primary;
assert isSortedBam secondary;

stage {
  name = "xenomapper-allocate";
  buildInputs = with pkgs; [ samtools python3Packages.xenomapper ];
  outputs = [ "out" ] ++ outs;
  buildCommand = ''
    xenomapper ${optionalString (flags != null) flags} \
      --primary_bam ${primary} --secondary_bam ${secondary} \
      ${concatMapStringsSep " " (out: "--${out} >(samtools view -bS - > ${"$" + out})") outs}

    mkdir -p $out
    for x in ${concatStringsSep " " outs} ; do
      ln -s ''${!x} $out/$x.bam
    done
  '';
}
