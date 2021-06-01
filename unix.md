# Introduction to Unix

This tutorial is based on one by Dr Nicholas Navin (MD Anderson Cancer Center, University of Texas), and another from Prof. Roger Butlin at University of Sheffield [here](https://openwetware.org/wiki/Butlin:Unix_for_Bioinformatics_-_basic_tutorial), as well as some additional material added by Dr Peri Bolton and Prof Chris Balakrishnan.


# Learning Objectives

The idea is for you to familiarise yourself with basic unix and some things we commonly use in bioinformatics. This will hopefully give you a solid foundation for the rest of this course, and for your future in genomic analysis. 
There are some questions embedded in the exercises. Please

* Understand and execute basic navigation of unix commandline
* Understand unix commandline terminology
* Understand and execute basic file and text manipulation
* Perform basic code troubleshooting
* Use aliases and variables
* write and execute shell scripts
* Ability to write basic loops

# Schedule

The first lesson will be part lecture (Wednesday 1100 EDT/0800 PDT), part interactive session that will orient you to the unix shell ([Directories &amp; some basic unix terminology.](#directories--some-basic-unix-terminology))
In your own time you will work through the remainder of the workshop, using the Slack Channel to ask questions. 
Then we will reconvene to go over the answers to the tutorial questions, as well as combat any additional questions you may have. (Thursday 1600 EDT/1300 PDT)


# Table of Contents

* [How To](#how-to)
	* [Things to remember:](#things-to-remember)
* [Introduction to Unix Tutorial](#introduction-to-unix-tutorial)
	* [Directories &amp; some basic unix terminology.](#directories--some-basic-unix-terminology)
	* [Jobs](#jobs)
	* [Doing things with files](#doing-things-with-files)
	* [Files and Fasta files](#files-and-fasta-files)
	* [File Manipulation](#file-manipulation)
* [Aliases](#aliases)
* [Variables](#variables)
* [Shell Scripts](#shell-scripts)
* [nohup](#nohup)
* [For Loops](#for-loops)


## How To

If you are on Windows, you need to open your ```WSL``` (e.g. Ubuntu)... Or if you have selected  ```PuTTY``` (please don't) you'll need to login to the server first. 

If you are on Mac you'll just need to open your ```Terminal```

### Things to remember:
* Pay attention to the helpful hints here and there.
* Google is also helpful (seriously, use it when you are stumped, I do it all the time) 
* For any command, you can always ask for information about the usage and options using the inbuilt manual (man) page from the commandline. E.g.,: ```man grep``` or ```grep --help```
* Unix is CasE SensiTive
* It is very easy to overwrite files if you name them the same thing. (You will not be asked “are you sure”?)
* Keep your directories organized or things spiral out of control quickly.
* A good text editor is nice to have (on macs, TextWrangler is great, and free). Unix geeks like Unix text editors nano and vi that can be run in Unix. 
* Even though I’ll show you how to transfer files using the command line, GUI ftp software are convenient (Cyberduck for macs, Filezilla, I think works on windows). (Again, a true Unix geek would frown on such things)
* you can use up and down arrows to scroll through previous commands 
* You can use ```tab``` to auto-complete names that are already in the system (This is useful for long filenames!)
* If you get stuck, Software Carpentry offers a neat [reference](https://swcarpentry.github.io/shell-novice/reference) for some basic unix commands. 

## Introduction to Unix Tutorial

Open your Mac Terminal or Linux WSL on Windows.

We are going to log into the university server using ```ssh```, which allows you to log into a remote computer. Please refer to your emails for the IP address and password. 

```bash
ssh ngsclass@<IP.ADDRESS>
```

You will be prompted for a password, type it in and click enter.


### Directories & some basic unix terminology.

#### 1.1.1 pwd


To determine the directory you are in use ```pwd``` which stands for **Print Working Directory**

You are logged into the home directory. This can be accessed by the full name ```/home/ngsclass/``` or using ```~/``` or ```$HOME```

#### 1.1.2 ls

You can **list** the current directory contents using ```ls``` or ```ls .``` where the ```.``` indicates the current directory. There will be many files in the home directory.

There is also a root directory. To view what is in that directory you can use ```ls /``` 

The ```ls``` **command** also has a number of **arguments** or **flags** that you can apply to it. For example ```ls -l```. Run it. Flags in unix are typically preceded by a ```-```

The output from a command is called **Standard Out** or ```stdout```, and is usually output on your screen. In-contrast, the information you give to a program to print or modify is called `stdin` (e.g. the `.` in `ls .`. If you make a mistake you will get an error, referred to as `stderr`.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

#### 1.1.3 mkdir & cd 

Now, make a directory for you for the duration of the course using ```mkdir```. Ideally this will be your last name, and should not match any folders already existing in the home directory. I have made one for myself and called it ```Bolton/```

For example:

```bash
mkdir Xena
```

now move into the directory you just made (note: don't just copy Xena)

```bash
cd Xena
ls
```

You've made an empty directory! 

**Tip:** You can use ```cd``` to specify a full directory path that is not directly connected to your current directory (e.g. ```cd /usr/local/bin/``` is where most programs are installed).

You can also access the previous directory you were in using ```cd -```

You can jump back a directory using ```cd ../``` or two ```cd ../../``` and so on.


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

#### 1.1.4 help pages

Most programs in unix have help pages associated with them. These are useful for knowing what flags are available to use, and what they do.
These are called ```man``` pages, meaning **manual** page. 

They can be invoked by invoking ```man``` plus the name of the program you want help with. For example:

```bash
man ls
```


**Question 1:** How do you get the file size to display in a more human readable way when using ```ls```? 


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### Jobs 

Let’s run the sleep program on the server.  It will put the computer to sleep for a designated number of seconds

```bash
sleep 10
```

As you have noticed this runs the process in the “foreground” which can prevent you from continuing to issue commands.  Instead, use the ‘&’ sign at the end of the line to run the process in the “background”.


```bash
sleep 260 &
```

If you start a job in the **foreground**  (without using ```&```) and want to stop it you can kill it completely using ```ctrl+c``` or suspend it using ```ctrl+z```.

Let's see what jobs your shell (session) has started. 

```bash
sleep 30
ctrl+z

jobs
```
	
You should be able to see two jobs - the ```[2] sleep 30``` and ```[1] sleep 260 &``` 
You can resume  ``` sleep 30``` and move it to the **background** by using ```bg %2```

We can see what processes are running using ```ps``` (**process status**)

You should still see ```sleep 260``` on the list of active processes. Use the Process ID (PID) in the first column to kill the process using the command ```kill <PID>```. Replacing everything including the brackets with the PID for the right proces.


Another useful command is `htop`. Not all servers have this installed, but it can be.  Run the command and see what it looks like.

`htop` is a live updated and interactive system-monitor and process manager. You can see what other users are doing on there, and how many cores are available, and how much RAM. This is useful for when the server is busy, and you can also use it to extract your PIDs and kill jobs. [Here's](https://linuxtogether.org/htop-command-explanation/) a basic intro to what `htop` shows.

**Question 2:** How many CPUs does this server have? How much total RAM is there? 



Other servers will have a job submission system such as `slurm` or `qsub`, which puts jobs into a queue and manages the server resources for the user. We will not get into these here. 


Now, let's play with some genome data.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### Doing things with files

First, make sure you're in the directory with your name on it that we made at the beginning of the tutorial. 

you can make files using ```touch```

```bash
touch enzyme.txt
```
you can also make multiple files at once

```bash
 touch test test1 test12 test123 test.file test.file.txt
```

Note that you cannot make a filename with a space in it (it will think it is two separate files). Safe characters for use in filenames are ```a-zA-Z0-9._```

Now create a directory called ```bacteria```. 

**Move** the ```enzyme.txt``` file into this folder using ```mv```

```bash
mv enzyme.txt
cd bacteria
ls
```
you could also **copy** the file into the directory using ```cp```

Now you can **remove** the file using ```rm```.  Please be careful with this command **There is no recycle bin**. Removed files cannot be restored. 

```bash
rm enzyme.txt
```

**Question 3:** Look at ```man rm``` page for arguments you can pass to rm. Why might you want to run ```rm -i <FILE>```?

now, move back to the parent directory for ```bacteria``` (the one with your name on it). 

you can delete this now empty directory with

```bash
rmdir enzyme.txt
```

you may also use ```rm -r```for a full directory or ```rm -d``` removes only empty directories

we can also use **wildcard** ```*``` to search for files that we created

```bash
ls -l test*
```
now remove all test files using ```rm``` and the wildcard.

#### echo, file contents and redirects.

the ```echo``` command is useful for doing all sort of things - particularly in sanity checks later when you build loops.

```bash
echo foo
```

By default, ```echo``` prints to ```stdout```.  You can use the ```>``` operator to redirect the content to a file.

```bash
echo foo > test1.txt
echo bar > test2.txt

cat test1.txt test2.txt
```

Notice, that the ```cat``` command has printed the contents of ```test2.txt``` right after the contents of ```test1.txt```. In other words it has con```cat```enated the two files.

The ```>>``` operator can be used to append content to an existing file.

```bash
echo "I am some text" >> test1.txt
cat test1.txt
```

Note that if you already have content in a file and you try to use ```>``` you will lose all previous content. **Try it**

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### Files and Fasta Files

Let's do some more basic commands using some genome data. We can download a .zip file that has been pre-prepared using ```wget``` or ```curl```

```bash
curl <GET URL FROM YOUR INSTRUCTOR> -o Archive.zip -J -L -k
```
[/]:# I am using a different version of this from what was originally provided in the Texas Tutorial. I am hosting on my Dropbox but won't forever.  

You should see a ```.zip``` file in your directory now. You need to ```unzip``` it

```bash
unzip Archive.zip
ls -l
```

We have some genome related files in this folder, including a gzipped **fasta** file. Fasta files are often indicated with the extensions: ```.fa``` ```.fna``` ```.fasta```
```gzip``` is another compression algorithm like `zip`.  

```bash
ls -l
gunzip chr1.fna.gzip
```

now we can take a sneak peak of the file contents

```bash
head chr1.fna
````

This is a typical fasta file. Sequence information is contained in the header preceded by the `>`. You can have multiple sequences in the same fasta file.

You can also look at the whole file at once using ```less```

```bash
less chr1.fna
```

you can scroll through the file using your `up` and `down` arrows. To exit you can either press `q` or `ctrl+z`. 

**Question 4:** What does the -S argument do? When might that be useful to you?

You can preview a file using `head` and `tail`.

```bash
head -n 25 chr1.fna
```

**Question 5:** What does the -n argument do?

You can also ask how many lines there are in a file.

```bash
wc -l chr1.fna
```

We can also **search** files for key phrases. For example: 
```bash
grep "chromosome" chr1.fna
```
We can also look at how many instances of a phrase are in a file using a **pipe**, which connects two commands together. Here we will look for the number of occurrences of the EcoRI restriction site (GAATTC).

```bash
grep "GAATTC" chr1.fna | wc -l
```

The pipe `|` takes the **Standard output** or `stdout` from the previous commands, and performs another command on it. In this case, we ran the search, then we used `wc -l` to count the number of instances of our search term.
We use `-l` because we know that `grep` prints the line that contains our search term. 

**Question 6:** Can you think of a reason why this is NOT the best way to find a restriction enzyme site (hint: think about what might happen at the end of line)?

We can also specify that we want to return everything but the search string. For example, we could count the number of nucleotides in the file.

```bash
grep -v ">" chr1.fna | wc -c 
```

**Question 7:** There is a file called `human_pep.fa`. How many protein sequences are there in this file? How many amino acids are there? Use a redirect to make a file with the names of all the proteins in this file.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### File Manipulation

#### Find and Replace

We can do 'find and replace' in a file. For example, here we want to **translate** a cDNA sequence into an RNA sequence (replace Thymine with Uracil) using `tr`. 

```bash
cat ZF_SNCA_cDNA.fa | tr "T" "U"
```


Another powerful replacement tool in UNIX is the sed (stream editor) command. As you can imagine, find/replace comes up A LOT. This command is useful. Let’s use sed to replace all the U’s with T’s in the `ZF_SNCA_cDNA.fa` file.
The format s/pattern1/pattern2/g is used, in which s indicates replace, and g stands for global.

```bash
sed 's/T/U/g' ZF_SNCA_cDNA.fa
```
**Question 8:** What is the problem with doing it these ways? How could you avoid this problem? I don't have the answer - time to get googling and see what works!

The `sed` command is more powerful than `tr` because it uses *‘regular expressions’*, patterns to match complex patterns and replace text, while the `tr` command is limited to replacing exact characters and words. With `sed` we can also delete patterns in text. Let’s delete all of the A’s in the kras fasta file.

```bash
sed 's/A//g' ZF_SNCA_cDNA.fa
```

**Question 9:** How would you use a *regular expression* pattern to replace all lines starting in A with an underscore character in the 'ZF_SNCA_cDNA.fa' file? Don't be afraid to google!

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

#### Splitting and Joining files

UNIX has some very useful tools for splitting large files into smaller ones. Let’s imagine that we want to split the first ten proteins in the human proteome file into individual FASTA files. We can use the csplit tool to accomplish this. We need to specify the name of each file, which we will call ‘file’ and the total number of sequences that we will split which is 10

```bash
csplit -k -f file human.pep.fa '/>/' {10}
ls
```

Now lets grab file numbers 1, 7 and 8 and merge them into a single file using `cat`

```bash
cat file01 file07 file08 > bigfile.txt
```

Now let’s check to make sure we did in fact merge three FASTA files

```bash
grep ‘>’ bigfile.txt | wc -l
```
<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

#### some manipulation of .bed files

BED files are used to store annotations or data values in a simple text format. BED files are related to numerous other file types that are common in bioinformatics/genomics/etc (gff, gff3, SAM). Namely, text columns are used to identify starting and stopping points of various stuctures (e.g., genes, introns) relative some reference (e.g., genome). BED files require at least three columns, with an optional column as a descriptor or gene name.

```
column1: chromosome number
column2: start position
column3: stop position
column4: gene name or identifier
```

Note: columns in a BED file are always tab-delimited

Let’s examine the `cancer_gene.bed` file with `head`

```bash
head cancer_genes.bed
```

Often it will be necessary to extract a subset of columns from a BED file to produce another file. Use the cut command to extract the gene identifier column(4) and make a new file with the gene names.

```bash	
cut -f 4 cancer_genes.bed > genes.txt
```

**Question 10:** Think back to our lesson yesterday, how might you achieve this in R?

**On Your Own** Examine the first lines of the `genes.txt`, to confirm that the 4th column was extracted.

The genes are out of order, let’s **sort** them alphabetically using the `sort` command

```bash
sort genes.txt
```

**Question 11:** Sort them in reverse alphabetical order. You might want to use `man` or google, to discover the various options for `sort`. How do you sort in reverse order?


Now let’s find all the genes in that contain the string ‘RAS‘
```bash
grep RAS genes.txt
```

**Question 12:** How many genes contain the string RAS ?


Now let’s return to the `cancer_genes.bed` file and use the `head` and `cut` commands to extract the first 10 lines from column 1 and output a file called `column1.txt`

```bash
head -10 cancer_genes.bed | cut -f 1 > column1.txt
```

Now let’s use head and cut to extract the first 10 lines from column 3 and output a file called column3.txt.

```bash
head -10 cancer_genes.bed | cut -f 3 > column3.txt
```

The counterpart of the `cut` command is the `paste` command which can be used to paste columns together. 

**On your own** Use the `paste` command to stitch together columns 1 and 3 and make a new file called `join_columns.txt`

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Aliases
some content from [here](https://shapeshed.com/unix-alias/)
So, you think typing `ls -l` all the time is annoying? Think that you can't be trusted with a naked `rm`?

You can create a **shell alias** to avoid typing out full commands all the time. 
To create an alias for the session you are in:

```bash
alias rm='rm -i'
``` 
**On Your Own** now try creating and removing a file. What happens?

However, when you close the terminal window you will lose this alias. 

You can use the `~/.bashrc` file to set an alias that will load every time. **LOOK (DO NOT MODIFY)** at what we have in the `~/.bashrc` for the server.

If you want to, you can set your own aliases.

### For Windows SL

Open a new WSL terminal window. You will be in your home directory. 

use `pico` to edit your `.bashrc` file.

Type in the alias you want to set (e.g. `alias ll='ls -lhF'`) and hit `enter` for every new alias you want to set.

Exit and save your file. 

Now, to update your current shell with the new alias run ```source .bashrc```. Everytime you load the shell you'll now have this new versions of `ls` available to you.

### For Mac

Open a new Terminal window or a new tab (just like a in a browser). You will be in your home directory. 

use `pico` to edit your `.bash_profile` file.

Type in the alias you want to set (e.g. `alias ll='ls -lhF'`) and hit `enter` for every new alias you want to set.

Exit and save your file. 

Now, to update your current shell with the new alias run ```source .bash_profile```. Everytime you load the Terminal you'll now have this new versions of `ls` available to you.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


## Variables

Sometimes we might want to set a **Variable**, which could be any character string, such as a filename, a random string of letters and numbers. These are held in the session and can be called upon later. Just like in R!

Usually, variables in the shell are written in uppercase, but they don't have to be. However, they can ONLY contain `a-zA-Z0-9_`. They also cannot be preceeded by a number. 

e.g. `VAR_1` is a valid variable name but `2_VAR` is not.

Now we can set a variable

```bash
VAR="hello world"
```

The variable can be accessed using the `$` character.

```bash
echo $VAR
```

You could also set a variable as a directory. For example:

```bash
DIR=/usr/local/bin
cd $DIR
```

When you're setting a directory as a variable, make sure you don't end with `/`. This can be important if you are using variables as shorthand and then appending filenames.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Shell Scripts

You can run a series of commands at once using a **shell script**. This is basically a text file that contains a bunch of shell comands that will be executed at once.

We can make a shell script by typing in to a text editor like `pico` or `vi`, or by continuously appending to a file. We will demonstrate the latter.

```bash
echo '#!/bin/sh' > shell_script.sh
echo 'echo Hello World' >> shell_script.sh
```

The `.sh` file extension is convention to indicate the type of file that it is. just like you can make a text file without `.txt`, but that is not very useful to someone coming in with no memory of what the files are.
Now use `cat` to look at the contents of your script. 

**Question 13:** Why are there two usages of `echo` in the above commands?

You can run the script with.

```bash
sh shell_script.sh
```
This executes the file in a new shell process.

There is another way to launch a shell script. First you need to make it **executable**. 

```bash
chmod 755 shell_script.sh
./shell_script.sh
```

This is just saying 'file in this directory'. But by making it executable the computer reads the file and the first line `#!/bin/sh` and knows how to execute it. This file could also be another form of executable (like many of the programs we will use), or another kind of script file.
For example `perl` scripts (often with end with `.pl`) are specified with the shebang line `#!/usr/bin/perl` which tells the computer to execute the file using perl.

**Note:** Anything the first line of your shell script that follows `#` is a **comment**. This is so you can annotate your code so you know what you did later. Just like in some of the R scripts yesterday.

**Fun-fact:** You may also see shebang lines with `#!/bin/bash` (that could be exectuted with `sh`). This is a slightly different type of unix shell, with slightly different commands, but are largely the same for most basic bioinformatics work (I've found). None of this is important today, but will be if you do complex unix programming it might be good to know the different shells `dash`, `ash`, and `pdksh`.
Usually if you use `sh` the Linux operating system should point to its preferred shell (e.g. Ubuntu `sh` points to `dash`). 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Nohup

Another cool trick is to set a command so that it will execute even if the terminal window (shell) isn't open anymore. For this, you can use `nohup`

Many processes in bioinformatics can take hours or days (or even weeks). So it's good to be able to close the window. You can run a bioinformatics command on the commandline

For example

```bash
nohup bwa mem ref.fa reads.fq > aln-se.sam &
```

Or you can run a series of commands that you have put within a shell script. 
**HOT TIP:** If you have a shell script where one command needs to be completed prior to executing another command, individual lines should be followed by `;`

For example:

```bash
#!/bin/sh
bwa index ref.fa ;
bwa mem ref.fa reads.fq > aln-se.sam ;
samtools view -bS aln-se.sam > aln-se.bam ;
```
Because these commands are sequential - you cannot run one without having the outputs from the previous file. 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## For-loops

You should remember the basic idea of for loops from yesterday's tutorial on control structures in R.
Basically, for *each* iteam in a list, do a thing.

Using a text editor like `pico` make the following for loop into a shell script. 

```bash
#!/bin/sh

ROTJ="Luke Leia Han Chewy R2 Ewoks Vader"

for j in $ROTJ
do
	echo $j
done
```

Name it something and run it. 

**Hot Tip:** You can specify a lot of files in a variable like above ... e.g. `FASTA=*.fa`


Another loop could be over a range of numbers

```bash
#!/bin/sh

for value in {1..10}
do
	echo $value
done
```

Run it!

**Question 14:** How would you write a script that counts the number of sequences in each of the fasta files we downloaded today?

Phew! That was a lot! 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>
