{ bionix
, flags ? null
} :

{ input1
, input2 ? null
} :

with bionix;
with lib;
with types;

let
  fq = f: matchFiletype "fastp-input" { fq = _: f; gz = matchFiletype' "fastp-input" { fq = _: f; }; } f;

  out =
    stage {
        name = "fastp";
        buildInputs = [ pkgs.fastp ];
        outputs = [ "out" "fastq1" "json" ] ++ (if input2 != null then [ "fastq2" ] else []);
        buildCommand = ''
            mkdir -p $out
            fastp \
                ${optionalString (flags != null) flags} \
                -i ${fq input1} \
                -o fastq1.fq.gz \
                ${optionalString (input2 != null) ''
                    -I ${fq input2} \
                    -O fastq2.fq.gz \

                    cp fastq2.fq.gz $fastq2
                ''}

                cp fastq1.fq.gz $fastq1
                cp fastp.html $out
                cp fastp.json $json
        '';
    };
  fqgz = { filetype = filetype.gz (filetype.fastq {}); };
in out // { fastq1 = out.fastq1 // fqgz; } // (if input2 != null then {fastq2 = out.fastq2 // fqgz; } else {})
