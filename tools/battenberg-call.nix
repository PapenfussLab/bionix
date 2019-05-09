{ bionix
  , gc-correction-loc
  , prob-loci
  , impute-info
  , thousand-genomes-loc
  , ignore-contigs-file
  , gender
  , assembly
  , species
  , flags ? null}:
  { normal, tumour }:
    stage {
      name = "battenberg";
      buildInputs = [
        (bionix.battenberg.app)
        ];
      buildCommand = ''
        mkdir $out
        battenberg.pl -o $out \
          -prob-loci ${prob-loci} \
          -impute-info ${impute-info} \
          -thousand-genomes-loc ${thousand-genomes-loc} \
          -ignore-contigs-file ${ignore-contigs-file} \
          -gender ${gender} \
          -gc-correction-loc ${gc-correction-loc} \
          -assembly ${assembly}\
          -species ${species} \
          -t $NIX_BUILD_CORES \
          ${lib.optionalString (flags != null) flags}
        '';
      passthru.multicore = true;
      }
