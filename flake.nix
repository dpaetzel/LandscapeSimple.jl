{
  inputs = {
    nixpkgs.url =
      # "github:dpaetzel/nixpkgs/dpaetzel/nixos-config";
      "github:dpaetzel/nixpkgs/update-clipmenu";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
    devenv.inputs.nixpkgs.follows = "nixpkgs";

    overlays.url = "github:dpaetzel/overlays/master";
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs =
    {
      self,
      nixpkgs,
      devenv,
      systems,
      overlays,
      ...
    }@inputs:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      devShells = forEachSystem (
        system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [
              overlays.overlays.mydefaults
            ];
          };
        in
        rec {
          default = devenv.lib.mkShell {
            inherit inputs pkgs;
            modules = [
              {
                # https://devenv.sh/reference/options/
                languages.julia.enable = true;
                languages.julia.package = pkgs.myjulia;

                scripts.repl.exec = ''
                  julia --project=.
                '';
              }
            ];
          };
          devShell.${system} = default;
        }
      );
    };
}
