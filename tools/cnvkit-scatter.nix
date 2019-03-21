{bionix
,flags ? null}:

input:

with bionix;
with lib;
with types;

stage {
  name = "cnvkit-scatter";
  buildInputs = [ cnvkit.app ];
  buildCommand = ''
    cnvkit.py scatter -s ${input}/*.cn{s,r} -o plot.pdf
    cp plot.pdf $out
  '';
}
