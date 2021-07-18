# #!/bin/bash
# Tips: `> output`


#################### Basic Pcp and Mt ####################
#[commands separate by lines and {}. end with ; in one sentence]
[command line1]; [command line2]; [command line3]

#################### How to make the shellscript with options ####################


#################### Mathematical Calculation ####################

var2=2var=$((1+1*2**$var2-12/2)) # the order is ok

#################### Assign variables ####################
set variablesvar = $(command to come up an object) # method1 (default and recommended)
echo ${var := "BPSM"} # method2
var2 = ${command with $var} # you could use the var to create another var

# use variables
echo $var 
# Use the self-defined var
echo ${var}"is fecth" 

# {} set the scale of the var. you can also use the $var in the commends and some situation should be considered.
# clear variables
unset var # No $ should be used with the unset command when pointing to var
# get the length of string variable
${#string}  
## e.g.
var="str"
echo ${#var}
> 3

# extract substring from the string variable
${string:position:length} 
## e.g.
var="stella"
echo ${var:4}
> stel
echo ${var:4:2}
> st

# delete the shortest/first total match of string variable
# Tips: append character every time to match the pattern and once match, then delete this part and return the reset of the string variable
${string#matchPattern}  # from front of the string
${string%matchPattern}  # from back of the string
## e.g
var="bpsm.get.txt"
var1="b.bpsm.get.txt"
echo ${var#*.}
> get.txt
echo ${var1#?.}
> bpsm.get.txt
echo ${var%.*}
> bpsm.get

# delete the longest/last total match of string variable
# Tips: the same but the last one would be considered as matched and delete
${string##matchPattern}  # from front of the string
${string%%matchPattern}  # from back of the string
# e.g
var="bpsm.get.txt"
var1="b.bpsm.get.txt"
echo ${var##*.}
> txt
echo ${var1%%.*}
> b
echo ${var%%.*}
> bpsm.get

# Replace string in string variable
${var/pattern/replacement}  # replace the first and longest match
${var//pattern/replacement}  # replace all match
${var/#pattern/replacement} # match from the front of the string
${var/%pattern/replacement} # match from the end of the string
# Tips: Remember the var is not changed by this action and you should assign this to antoher variable.
# e.g.
var="go.exe.exe.zip.exe"
echo ${var/.e*./}
> goexe
echo ${var/.e*./.tar.}
> go.tar.exe
echo ${var//.e*./.tar.}
> go.tar.exe
var="go.exe.tarexe.zip.exe"
echo ${var//.exe./.tar.}
> go.tar.tarexe.zip.exe

#################### Control Structure ####################
########## Judge ##########
## if 
if [condition]
then
	[commands]
elif
	[commands]
else	
    [commands]
fi # it is necessary for the bash shell

########## Loop ##########
# for loop
for x in [sequences]; do[commands]  # temp var would be specified as $x done
## e.g.
for x in {a..z}; 
do 
egrep -o 'D. '$x'\S+' dros_list.txt > $x'_species.txt' 

# while loop
## Tips: the variables are separated by the IFS
while read [local variables] 
do
[contents]
done < [input]
## e.g.
1ls -1 *.txt > myfiles.list
cat myfiles.list
while read myfilename
do
echo -e "Processing $myfilename..."
linesinfile=$(wc -l $myfilename | cut -d ' ' -f1)
echo -e "\thas $linesinfile lines in it"
done <  myfiles.list#
# e.g.2 count the line in different documents
while read myfilename
do
echo -e "Processing $myfilename..."
linesinfile=$(wc -l $myfilename | cut -d ' ' -f1)
echo -e "\thas $linesinfile lines in it"
done <  myfiles.list

#################### Fucntions ####################
function x(){	
	# $1 and $2 will follow the input order
	$1 + $2
	}# 


#################### AWK Methods ####################
#!/usr/bin/awk -f
BEGIN {FS="\t";OFS="_";} {count++;if($1 != "#")  
{
	print "Currently doing "count ;
	total=total + ($12 * $3)  }}
	END {print "The total for "count" lines was " int(total) > "awkoutputfile.txt";
	print "Script run complete." >> "awkoutputfile.txt" ;
	print "Script run complete."system("ls -alrt *awk*")
}

