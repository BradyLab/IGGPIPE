# This converts a .asciidoc file into a .html file.
#   $1: Name of file to convert, without .asciidoc suffix that file must have
#
# This requires that:
#   1. "python" is installed on your system (https://www.python.org/downloads/)
#   2. "asciidoc" tools are installed (http://asciidoc.org/INSTALL.html)
# Make a symbolic link in a directory on your path that is named asciidoc,
# pointing to the "asciidoc.py" file in the asciidoc install directory.
# For example, if asciidoc-8.6.9 were downloaded into ~/src, and if ~/bin
# is on your path, this would make a symbolic link:
#   ln -s ~/src/asciidoc-8.6.9/asciidoc.py ~/bin/asciidoc
# Also make sure the file "asciidoc.py" is executable:
#   chmod +x ~/src/asciidoc-8.6.9/asciidoc.py
#
asciidoc -b html $1.asciidoc
