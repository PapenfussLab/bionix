{ stdenv
, callPackage
, lib
, bwa
, samtools ? null
, ref
, bamOutput ? true
, flags ? null
}:

{ input1
, input2 ? null
}:

assert bamOutput -> samtools != null;

with lib;

let index = callPackage ./bwa-index.nix { inherit bwa stdenv lib; } ref;

in stdenv.mkDerivation {
  name = "bwa-mem";
  buildInputs = [ bwa ] ++ optional bamOutput samtools;
  buildCommand = ''
    ln -s ${ref} ref.fa
    for f in ${index}/* ; do
      ln -s $f
    done
    bwa mem ${optionalString (flags != null) flags} -t $NIX_BUILD_CORES ref.fa ${input1} \
      ${optionalString (input2 != null) input2} \
      ${optionalString bamOutput "| samtools view -b"} \
      > $out
  '';
}
