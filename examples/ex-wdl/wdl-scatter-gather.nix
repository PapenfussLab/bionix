# The scatter-gather example from https://github.com/openwdl/wdl
# translated to bionix
{ bionix ? import ./../.. { } }:

with bionix;
with lib;

let

  prepare = splitString "\n" (removeSuffix "\n" (readFile (stage {
    name = "prepare";
    buildInputs = [ pkgs.python3 ];
    buildCommand = ''
      python -c "print('one\ntwo\nthree\nfour', end=''')" > $out
    '';
  })));

  analysis = str: removeSuffix "\n" (readFile (stage {
    name = "analysis";
    buildInputs = [ pkgs.python ];
    buildCommand = ''
      python -c "print('_${str}_')" > $out
    '';
  }));

  gather = strs: stage {
    name = "gather";
    buildCommand = ''
      echo ${concatStringsSep " " strs} > $out
    '';
  };

in
gather (map analysis prepare)
