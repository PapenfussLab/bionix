{
  description = "Bioinformatics workflows with Nix";
  outputs = { self, ... }: {
    lib = import ./.;
  };
}
