hisat2 example:
hisat2-build $hisat2/example/reference/22_20-21M.fa --snp $hisat2/example/reference/22_20-21M.snp 22_20-21M_snp
hisat2 -f -x $hisat2/example/index/22_20-21M_snp -U $hisat2/example/reads/reads_1.fa -S eg1.sam
head eg1.sam
hisat2 -f -x $hisat2/example/index/22_20-21M_snp -1 $hisat2/example/reads/reads_1.fa -2 $hisat2/example/reads/reads_2.fa -S eg2.sam
samtools view -bS eg2.sam > eg2.bam
samtools sort eg2.bam -o eg2.sorted.bam
samtools mpileup -uf $hisat2/example/reference/22_20-21M.fa eg2.sorted.bam | bcftools call -c -v --output-type b  | bcftools view -m 2 --output-type b - > eg2.raw.bcf
bcftools view eg2.raw.bcf




name a path with hisat2 and check it
export hisat2=/project/bioinformatics/Kim_lab/s179814/hisat2
echo $hisat2



in emacs python mass-comment/uncomment 
C-x r t cursor should be in the first column of the last row
C-x r k cursor should be in the second column of the last row


>wc -l <filename>
count the number of lines in a text file

>sed -n '4149273,4149300p;4149300q' genome.fa  > genome.cut
copy lines from genome.fa and insert in genome.cut

read text file from the first few lines 
less <filename>

Bash edit
emacs -nw ~/.bashrc

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

C-s	Search
M-%	Replace

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
git add -A
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

