{
  bionix,
  gurobiLicense ? null,
  count_reads ? true,
  genotype_snps ? true,
  count_alleles ? true,
  combine_counts ? true,
  cluster_bins ? true,
  plot_bins ? true,
  compute_cn ? true,
  plot_cn ? true,
  size ? "50kb",
  mincov ? 8,
  maxcov ? 300,
  snps,
  phase ? "None",
  diploidbaf ? 0.08,
  tolerancerdr ? 0.15,
  tolerancebaf ? 0.04,
  sizethreshold ? 0.01,
  figsize ? "6,3",
  clones ? [2 6],
  seeds ? 400,
  minprop ? 0.03,
  diploidcmax ? 6,
  tetraploidcmax ? 12,
  ghostprop ? 0.35,
  limitinc ? 0.6,
  blocklength ? "50kb",
}:
with bionix;
with lib;
with types;
  {
    normal,
    tumours,
  }: let
    getRef = matchFiletype "hatchet" {bam = {ref, ...}: ref;};
  in
    assert all (x: getRef normal == getRef x) tumours; let
      ref = getRef normal;

      ini = pkgs.writeText "hatchet.ini" (generators.toINI {} {
        run = {
          inherit count_reads genotype_snps count_alleles combine_counts cluster_bins plot_bins compute_cn plot_cn;
          reference = "${lnRef ref}/ref.fa";
          normal = "${lnBam normal}/input.bam";
          bams = concatMapStringsSep " " (x: "${lnBam x}/input.bam") tumours;
          samples = "@SAMPLES@";
          output = "./out";
          processes = "@PROCESSES@";
        };

        count_reads = {inherit size;};
        genotype_snps = {
          inherit mincov maxcov;
          snps = "${lnVcfBz snps}/vcf.bgz";
        };

        count_alleles = {inherit mincov maxcov;};
        combine_counts = {inherit blocklength phase;};
        cluster_bins = {inherit diploidbaf tolerancerdr tolerancebaf;};
        plot_bins = {inherit sizethreshold figsize;};
        compute_cn = {
          inherit seeds minprop diploidcmax tetraploidcmax ghostprop limitinc;
          clones = concatMapStringsSep "," builtins.toString clones;
        };
      });

      lnRef = ref:
        linkOutputs {
          "ref.fa" = ref;
          "ref.fa.fai" = samtools.faidx {} ref;
          "ref.dict" = samtools.dict {} ref;
        };

      lnBam = bam:
        linkOutputs {
          "input.bam" = bam;
          "input.bam.bai" = samtools.index {} bam;
        };

      lnVcfBz = vcf: let
        bz = compression.bgzip {} vcf;
      in
        linkOutputs {
          "vcf.bgz" = bz;
          "vcf.bgz.tbi" = samtools.tabix {} bz;
        };

      getSN = let
        script = pkgs.writeText "getSN.awk" ''
          BEGIN{
            FS=":"
            RS="[ \t\n]"
          }
          $1=="SM"{print $2; exit}
        '';
      in
        pkgs.writeShellScriptBin "getSN" ''
          exec samtools view -H "$1" | awk -f ${script}
        '';
    in
      stage
      ({
          name = "HATCHet";
          buildInputs =
            [hatchet.app getSN pkgs.samtools pkgs.bcftools]
            ++ (
              if gurobiLicense == null
              then [pkgs.cbc]
              else [pkgs.gurobi]
            );
          buildCommand = ''
            # Get tumour names
            names="${concatMapStringsSep " " (x: "$(getSN ${x})") tumours}"

            # macro substitute ini file
            substitute ${ini} hatchet.ini \
              --replace "@PROCESSES@" "$NIX_BUILD_CORES" \
              --replace "@SAMPLES@" "$names"

            hatchet run hatchet.ini

            cp -r out $out
          '';
          passthru.multicore = true;
        }
        // (
          if gurobiLicense != null
          then {
            GRB_LICENSE_FILE = gurobiLicense;
            HATCHET_COMPUTE_CN_SOLVER = "gurobi";
          }
          else {
            HATCHET_COMPUTE_CN_SOLVER = "cbc";
          }
        ))
