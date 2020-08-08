cp main_8g.cxx main.cxx
make 
./a.out &> .output
sed -e '/#/d'  .output  > cbz.txt  
rm .output