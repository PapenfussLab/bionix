{ bionix
, indexFlags ? { }
, bias ? false
, bootstrapSamples ? 0
, seed ? 42
, plaintext ? false
, fusion ? false
, single ? false
, frStranded ? false
, rfStranded ? false
, fragmentLength ? null
, fragmentSD ? null
, ref
}:

with bionix;
with lib;

assert (!single || (fragmentLength != null && fragmentSD != null));

inputs:

let
  inherit (bionix.types) matchFiletype';
  isFastQ = matchFiletype' "kallisto-quant" { fq = _: true; gz = isFastQ; };

  empty = ./kallisto-quant-empty.h5;

  python = pkgs.python3Packages.python.withPackages (p: with p; [ h5py ]);

  noStamp = pkgs.writeScript "nostamp.py" ''
    #!${python}/bin/python
    import h5py
    import os
    def copy(obj, out, path):
      if type(obj) in [h5py._hl.group.Group,h5py._hl.files.File]:
        for key in obj.keys():
          if key != "start_time":
            copy(obj[key], out, path + "/" + key)
      elif type(obj)==h5py._hl.dataset.Dataset:
        out.create_dataset(path, data=f[path], track_order=False, track_times = False)
    with h5py.File(os.environ['out'] + "/abundance.h5", "r") as f:
      with h5py.File("repack.h5", "a", track_order = False) as g:
        copy(f, g, "")
  '';

in

assert (all (x: isFastQ x.filetype) inputs);

stage {
  name = "kallisto-quant";
  buildInputs = with pkgs; [ kallisto ];
  buildCommand = ''
    mkdir $out
    kallisto quant \
      -i ${bionix.kallisto.index indexFlags ref} \
      -o $out \
      ${optionalString bias "--bias"} \
      ${optionalString (bootstrapSamples > 0) "-b ${toString bootstrapSamples} --seed=${toString seed}"} \
      ${optionalString plaintext "--plaintext"} \
      ${optionalString fusion "--fusion"} \
      ${optionalString single "--single -l ${toString fragmentLength} -s ${toString fragmentSD}"} \
      ${optionalString frStranded "--fr-stranded"} \
      ${optionalString rfStranded "--rf-stranded"} \
      -t $NIX_BUILD_CORES \
      ${concatStringsSep " " inputs}

    # Make deterministic by removing timestamps and using hdf5 empty template
    cp ${empty} repack.h5
    chmod 644 repack.h5
    ${noStamp}
    cp repack.h5 $out/abundance.h5
    sed -i $out/run_info.json -e '/"start_time"/d'
    sed -i $out/run_info.json -e '/"call"/d'
  '';
  passthru.multicore = true;
}
