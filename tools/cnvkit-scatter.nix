{bionix
,flags ? null}:

input:

with bionix;
with lib;
with types;

stage {
  name = "cnvkit-scatter";
  buildInputs = with pkgs; [ pkgs.python3Packages.cnvkit ];
  buildCommand = ''
    cnvkit.py scatter -s ${input}/*.cn{s,r} -o plot.pdf
    cp plot.pdf $out
  '';
}
