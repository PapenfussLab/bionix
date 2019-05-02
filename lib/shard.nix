{bionix}:

with bionix;
with lib.types;

{
  regex = n: (exec (_: input: let
    re = let f = matchFiletype "shard-regex" {
        fa = _: "^>";
        fq = _: "^@";
        gz = f;
        bz2 = f;
      };
    in f input;
    decompress = matchFiletype "shard-regex-decompression" {
        fa = _: "cat";
        fq = _: "cat";
        gz = _: "gunzip";
        bz2 = _: "bunzip2";
      } input;
    compress = matchFiletype "shard-regex-compression" {
        fa = _: "cat";
        fq = _: "cat";
        gz = _: "gzip";
        bz2 = _: "bzip2";
      } input;
    compressPkgs = with bionix.pkgs; matchFiletype "shard-regex-compression" {
        fa = _: [];
        fq = _: [];
        gz = _: [ gzip ];
        bz2 = _: [ bzip2 ];
      } input;
  in outputDrvs (stage {
    name = "shard";
    outputs = [ "out" ] ++ builtins.genList (i: "out" + toString (i + 2)) (n - 1);
    buildInputs = [ pkgs.gawk ] ++ compressPkgs;
    buildCommand = let
      awkScript = pkgs.writeText "shard.awk" ''
        BEGIN{cout=0}
        FNR==NR{out[nout++] = $0;next}
        /${re}/{cout = (cout + 1) % nout}
        {print > out[cout]}
      '';
    in ''
      for o in $outputs ; do
        echo $(basename ''${!o}) >> outputs
      done
      awk -f ${awkScript} outputs <(${decompress} < ${input})
      for o in $outputs ; do
        ${compress} < $(basename ''${!o}) > ''${!o}
      done
    '';
    passthru.filetype = input.filetype;
    }))) [];

  fastQPair = n: {input1, input2}: lib.zipListsWith (i: j: {input1 = i; input2 = j;}) (lib.shard.regex n input1) (lib.shard.regex n input2);
}
