{ bionix
, flags ? null
} :

{ input1
, input2 ? null
} :

with bionix;
with pkgs.lib;

# Match input file type—how to do .fq and .fq.gz? Does bz2 work?

stage {
    name = "fastp";
    buildInputs = [ fastp.app ];
    outputs = [ "out" "fastq1" "fastq2" "html" "json" ];
    buildCommand = ''
        mkdir -p $out
        fastp \
            ${optionalString (flags != null) flags} \
            -i ${input1} \
            -o fastq1.fq.gz \
            ${optionalString (input2 != null) ''
                -I ${input2} \
                -O fastq2.fq.gz \

                cp fastq2.fq.gz $fastq2
                ln -s $fastq2 $out/fastq2.fq.gz
            ''}

            cp fastq1.fq.gz $fastq1
            cp fastp.html $html
            cp fastp.json $json

            ln -s $fastq1 $out/fastq1.fq.gz
            ln -s $html $out/fastp.html
            ln -s $json $out/fastp.json
    '';
}