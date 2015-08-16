# Courtesy of user JPS at apple.stackexchange.com, 17-Jul-2015.
application=$(basename "$0")
if [ "$#" == 0 ]; then
    echo "Usage: $application path [paths...]"
    exit 1
fi
trashdir="/Users/${USER}/.Trash"
while (( "$#" )); do
    if [ -e "$1" ]; then
        src=$(basename "$1")
        dst=$src
        while [ -e "$trashdir/$dst" ]; do
            dst=$src+`date +%H-%M-%S`
        done
        mv -f "$1" "$trashdir/$dst"
    else
        echo "$1" does not exist.
    fi
    shift
done
