#!/bin/bash
#SBATCH -J MY_JOB_NAME
#SBATCH -p small
#SBATCH --ntasks-per-node=4
#SBATCH -o Oerr/%x-%j.out
#SBATCH -e Oerr/%x-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2021-08-06 21:24:45
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-09-04 13:40:51
 * @FilePath: /metaSC/Bash/relapwd.sh
 * @Description:
!EOF!


path_to=`readlink -m $1`

if [[ -z "${2}" ]]; then
    path_from=`pwd`
else
    path_from=`readlink -m $2`
fi

path_to=(${path_to//\// })
path_from=(${path_from//\// })
for element in $(seq 0 $((${#path_from[@]} - 1)))
do
        if [[ ${path_from[$element]} == ${path_to[$element]} ]];then
                path_from[$element]=""
                path_to[$element]=""
        else
                break
        fi
done

path_from=($(echo ${path_from[@]}))

if [[ -n ${path_from} ]];then
        for element in $(seq 0 $((${#path_from[@]} - 1)))
        do
                path_from[$element]="../"
        done
else
        path_from[0]="./"
fi
ddd=$(echo ${path_from[@]})
kkk=$(echo ${path_to[@]})

new_path=${ddd// /}${kkk// //}
echo $new_path
