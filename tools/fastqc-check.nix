{ bionix
, flags ? null
}:

with bionix;
with lib;

input:

stage {
  name = "fastqc-check";
  buildInputs = [ bionix.fastqc.fastqc pkgs.unzip ];
  stripStorePaths = false; # we do it explicity for fastqc
  outputs = [ "out" "zip" ];
  buildCommand = ''
    fastqc \
      -o $TMPDIR \
      ${optionalString (flags != null) flags} \
      ${input}

    sed  "s|$(basename ${input})|input|g" *.html > $out
    cp *.zip $zip
  '';
}
