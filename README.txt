wget http://cbcb.umd.edu/confcour/CMSC828H-materials/lab1.tar.gz
tar xvzf lab1.tar.gz

To connect to a server computer
ssh nucleus

Emacs operations

To edit a file
emacs -nw <filename>

To save the current file
ctrl+x,ctrl+s

To exit Emacs
ctrl+x,ctrl+c

To cancel pending commands
ctrl+g

To create a python script,
emacs -nw lab1.py, implement, exit

To run a python script
chmod +x lab1.py
./lab1.py

To undo
C-/
C-x u
C-_



To view a text file in the terminal
cat <filename>
less <filename> space for next page, q for quit

To clone a Git repository
git clone https://github.com/leejj82/compbio

git add <filename>
git commit -m "."
git push

git pull

awk '{if($3 >= 100) print}' NC_002951.ptt | wc -l

module (add programs such as python 3)
module avail
module list
module add
module rm

python (check python version)

c++ compile

>make
>./lab1