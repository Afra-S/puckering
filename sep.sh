

gmx trjconv -f md_red.xtc -s start.pdb -skip 200 -sep -o sep_snap_.pdb<<EOF
0
EOF

for file in `ls sep_snap_*`
do
    awk 'BEGIN{FS="";OFS="";prev_res_num=0;chain="A"}{res_num=substr($0,23,4);if($0~/^ATOM/ && res_num<prev_res_num){chain="B"} if ($0~/^ATOM/){$22=chain};print $0;prev_res_num=res_num}'  $file > sep_snap_temp 
    awk 'BEGIN{FS=""}{if($22=="B"){print $0}}' sep_snap_temp  > lig$T.pdb
    echo $T lig$T.pdb
    T=`expr $T + 1 `
done > list_structures.temp
