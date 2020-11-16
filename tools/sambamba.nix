{ bionix }:

with bionix;

let
  gen = callBionixE ./sambamba-generic.nix;

in {
  sort = callBionixE ./sambamba-sort.nix;
  index = def gen {tool = "index"; };
  merge = def gen {tool = "merge"; };
  slice = def gen {tool = "slice"; };
  flagstat = def gen {tool = "flagstat"; };
  markdup = def gen {tool = "markdup"; };
}
