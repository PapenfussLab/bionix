{ bionix
, nixpkgs
, ref
, bamOutput ? true
, flags ? null
, indexAttrs ? {}
}:

{ input1
, input2 ? null
}:

with nixpkgs;
with lib;

stdenv.mkDerivation {
  name = "bwa-mem";
  buildInputs = [ bwa bc ] ++ optional bamOutput samtools;
  buildCommand = ''
    ln -s ${ref} ref.fa
    for f in ${bionix.bwa.index indexAttrs ref}/* ; do
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
