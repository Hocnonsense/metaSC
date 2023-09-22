#!/bin/bash
#SBATCH -J MY_JOB_NAME
#SBATCH -p small
#SBATCH --ntasks-per-node=4
#SBATCH -o log/%x-%j.out
#SBATCH -e log/%x-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2023-08-26 20:45:01
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-22 20:45:06
 * @FilePath: /metaSC/Bash/delink.sh
 * @Description:
    change the linked path to a real path, and keep file in by link again
    including both file and hide files
!EOF!
set -e

function show_usage (){
    echo ""
    echo "Usage: $(basename $0) [-a] path"
    echo "Required arguments:"
    echo "    path | path to change from link to folder"
    echo "Optional arguments:"
    echo "    -a | if set, hide file will be linked as well"
    echo ""
    echo "FILE output:"
    echo "   Before:"
    echo "      a -> b"
    echo "      b/"
    echo "          c"
    echo "   After:"
    echo "      a/"
    echo "          c -> /real path to/b/c"
    echo ""

    return 0
}
declare use_a=""
declare path=""


set -- `getopt ha "$@"`

while [ -n "$1" ]
do
    case "$1" in
    -h)
        show_usage
        exit ;;
    -a)
        use_a="-a"
        shift ;;
    --) ;;
    *)
        path=$1
        break ;;
    esac
    shift
done

if [ ! -e $path ]
then
    show_usage
    echo "No path given, abort"
    exit 1
fi


declare real_path=`realpath $path`

if [ ! -L "$path" ]
then
    echo "path is not a link, abort" >&2
    exit 1
fi
if [ -f "$real_path" ]
then
    echo "path is not a dir, skip" >&2
    exit 0
fi

rm "$path"
mkdir "$path"

echo update $path with $real_path >&2

for j in `ls "$real_path" $use_a`
do
    if [ -d $path/"$j" ]
    then
        echo skip "$j"
    elif [ -e $path/"$j" ]
    then
        ln -s "$real_path"/"$j" "$1"/"$j"
        echo find broken link "$real_path"/"$j"
    else
        ln -s "$real_path"/"$j" "$1"/"$j"
        touch -achmr `realpath "$real_path"/"$j"` "$1"/"$j" | printf ""
    fi
done
