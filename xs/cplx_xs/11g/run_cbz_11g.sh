cp main_11g.cxx main.cxx
make 
./a.out &> .output
sed -e '/#/d'  .output  > flux_cbz.txt  
rm .output