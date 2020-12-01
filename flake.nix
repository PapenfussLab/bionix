{
  description = "BioNix";

  outputs = { self }: rec {
    lib = {
      bionixConfiguration = { configuration
        , pkgs }:
        import ./modules.nix {
          inherit pkgs;
          configuration = { ... }: { imports = [ configuration ]; };
        };
    };
  };
}
