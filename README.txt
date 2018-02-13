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

To cut the text
C-w .
To copy the text
M-w .
To paste the text
C-y .

To move back to screen fg




mkdir leejj
git clone https://github.com/leejj82/compbio
cp ~/work/lab1.cpp
git diff
git config --global user.email "leejj82@gmail.com"
git config --global user.name "leejj82"
git commit -a -m "."
git push
git pull

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

>time ./lab01 > /dev/null

