# This converts a .asciidoc file into a .pdf file.
#   $1: Name of file to convert, without .asciidoc suffix that file must have
#
# This requires that:
#   1. "python" is installed on your system (https://www.python.org/downloads/)
#   2. "asciidoc" tools are installed (http://asciidoc.org/INSTALL.html)
#   3. "fop" is installed (https://xmlgraphics.apache.org/fop/download.html).
# Use a "source download" for fop.  Make symbolic links in a directory on your
# path that are named asciidoc, a2x, and fop, pointing to the "asciidoc.py"
# and "a2x.py" files in the asciidoc install directory and to the "fop" file
# in the fop install directory.  For example, if fop-2.0 were downloaded into
# ~/java, and if ~/bin is on your path, this would make a symbolic link:
#   ln -s ~/java/fop-2.0/fop ~/bin/fop
# Also make sure the files "asciidoc.py", "a2x.py", and "fop" are executable:
#   chmod +x ~/java/fop-2.0/fop
#
asciidoc -b docbook $1.asciidoc
a2x -f pdf --fop $1.xml
