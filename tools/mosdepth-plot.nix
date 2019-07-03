{ bionix
, what ? "global"
, flags ? null
}:

with bionix;
with lib;

{ inputs
, names ? []}:

stage {
  name = "mosdepth-depth";
  buildInputs = with pkgs; [ python ];
  buildCommand = ''
    python ${pkgs.mosdepth.src}/scripts/plot-dist.py ${concatMapStringsSep " " (i: "${i}/out.mosdepth.${what}.dist.txt") inputs}
    ${concatStringsSep "\n" (zipListsWith (i: n: "substituteInPlace dist.html --replace ${i}/out.mosdepth.${what}.dist.txt ${n}") inputs names)}
    substituteInPlace dist.html --replace "x: 0.1" "x: 0.9"
    substituteInPlace dist.html --replace "y: 0.1" "y: 0.9"
    mv dist.html $out
  '';
}
