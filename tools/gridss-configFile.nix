{bionix}:

with bionix;
with lib;

let
    attrsToGridssConfigString = attrsToGridssConfigStringPrepend "";

    attrsToGridssConfigStringPrepend = prepend: attrs:
        concatStringsSep "\n" (
            attrValues (
                mapAttrs
                    (name: attr: prepend + (iniLine name attr))
                    attrs));

    iniLine = name: attr:
        let attrType = builtins.typeOf attr;
        in
            if (iniLineByAttrType ? ${attrType})
            then (iniLineByAttrType.${attrType} name attr)
            else builtins.throw (
               "`gridssConfig` cannot convert attribute of type \"" + attrType + "\".");

    iniLineByAttrType = {
        string = name: attr: name + " = " + attr;
        int    = name: attr: name + " = " + builtins.toString attr;
        float  = name: attr: name + " = " + (
                    builtins.head (
                        builtins.match "([0-9]+\.0?[1-9]*)0+" (builtins.toString attr)));
        bool   = name: attr: name + " = " + (if attr == true then "true" else "false");
        set    = name: attr: attrsToGridssConfigStringPrepend (name + ".") attr;
        # Allows for repeated fields (e.g. for adapters):
        list   = name: attr: concatStringsSep "\n" (map (x: iniLine name x) attr);
    };
in configAttrs: (pkgs.writeText
        "gridss.properties.override"
        ((attrsToGridssConfigString configAttrs) + "\n"))
