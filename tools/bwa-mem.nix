{ stdenv
, callPackage
, lib
, bc
, bwa
, index ? callPackage ./bwa-index.nix { inherit bwa stdenv lib; }
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

stdenv.mkDerivation {
  name = "bwa-mem";
  buildInputs = [ bwa bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    ln -s ${ref} ref.fa
    for f in ${index ref}/* ; do
      ln -s $f
    done
    cores=$(echo $NIX_BUILD_CORES ${optionalString bamOutput "- 1"} | bc)
    if [[ $cores -lt 1 ]] ; then
      >&2 echo "not enough build cores"
      exit 1
    fi
    bwa mem ${optionalString (flags != null) flags} -t $cores ref.fa ${input1} \
      ${optionalString (input2 != null) input2} \
      ${optionalString bamOutput "| samtools view -b"} \
      > $out
  '';
}
