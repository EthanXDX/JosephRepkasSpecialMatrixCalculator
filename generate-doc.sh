pandoc --pdf-engine=lualatex -V fontsize=12pt -H head.tex --highlight-style tango -f commonmark README.md -o README.pdf
