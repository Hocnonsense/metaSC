declare path=$1
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
 * @LastEditTime: 2023-08-26 20:46:08
 * @FilePath: /2023_08-M_mem-release/workflow/utils/libs/metaSC/Bash/delink.sh
 * @Description:
    change the linked path to a real path, and keep file in by link again
!EOF!
set -e

declare real_path=`realpath $1`

if [[ $path==$real_path ]]
then
    echo "path is not a link, abort" >&2
    exit 1
fi
if [[ -f $real_path ]]
then
    echo "path is not a path, skip" >&2
    exit 0
fi

rm $1
mkdir $1
ln -s $real_path/* $1

for i in $i
do touch -achmr `realpath $i` $i
done
