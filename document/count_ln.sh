line_num=0
fbase="../eostool/"
for file in `ls ${fbase}*pp`
do
    if [ ${file} != "${fbase}data.hpp" ]
    then
        new_ln=`wc ${file} | awk '{print $1}'`
        line_num=`expr ${line_num} + ${new_ln}`
        echo "File: \033[4;36m ${file} \033[0m, line number: \033[5;47;34m ${new_ln} \033[0m, total line number: \033[1;41;33m ${line_num} \033[0m "
    fi
done