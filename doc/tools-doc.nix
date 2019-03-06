{ bionix ? import ./.. {} }:

with bionix;

stage {
  name = "tools-docs";
  src = ../tools;

  xsltFlags = lib.concatStringsSep " " [
    #"--param section.autolabel 1"
    #"--param section.label.includes.component.label 1"
    #"--stringparam html.stylesheet 'style.css overrides.css highlightjs/mono-blue.css'"
    #"--stringparam html.script './highlightjs/highlight.pack.js ./highlightjs/loader.js'"
    #"--param xref.with.number.and.title 1"
    #"--param toc.section.depth 3"
    #"--stringparam admon.style ''"
    #"--stringparam callout.graphics.extension .svg"
  ];


  buildInputs = with pkgs; [ nixdoc libxslt libxml2 ];
  installPhase = ''
    function docgen {
      nixdoc -c "$1" -d "$2" -f "$1.nix" | sed 's/lib\./bionix./g' |grep -v locations.xml > "$1.xml"
    }

    docgen ascat 'ascatNGS CNV caller'
    docgen bowtie 'Bowtie aligner'
    docgen bwa 'BWA aligner'
    docgen cnvkit 'CNVkit CNV caller'
    docgen facets 'Facets CNV caller'
    docgen strelka 'Strelka2 variant caller'

    mkdir $out
    cp ${./tools.xml} tools.xml
    xmllint --nonet --xinclude --noxincludenode tools.xml --output tools-full.xml
    cat tools-full.xml
    xsltproc $xsltFlags \
      --nonet \
      --xinclude \
      --output $out/index.html \
      ${pkgs.docbook_xsl_ns}/xml/xsl/docbook/epub/docbook.xsl \
      tools-full.xml
  '';
}
