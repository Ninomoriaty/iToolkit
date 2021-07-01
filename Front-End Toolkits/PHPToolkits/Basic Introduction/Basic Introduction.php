<html lang="en">
<body>
<pre>
<?php
# Search for Manual
echo "\n------------------ Help and Manual ------------------\n";



# Comments
echo "\n------------------ Comments ------------------\n";
// Comments for one line
/* Comments
for
multiline
*/


# Variables
echo "\n------------------ Variables ------------------\n";
## Assign value to variables
$var = "assign method";

# Data Types
## Associative Array
echo "\n--- Associative Array ---\n";


## String
echo "\n--- String ---\n";
$str0 = "string";
$str1 = 'string';
### String Manipulation
### Tips: Double quotes allow variable substitution within the quotes and single quotes do not.
echo $str0 . "\n";
echo $str0 . "$var" . "\n";
echo $str0 . '$var' . "\n";
#### Join strings with space-separated . symbol
echo "string 1 link to " . "string2" . "\n";
### Upload textarea
/*
<form action=“page2.php" method="post"><pre>
    Your ID  <input type="text" name="userid"/>
    <textarea cols=80, rows=25 name=“mytext"/>
    </textarea>
    <input type="submit" value="go" />
    </pre></form>

And handle with the following PHP
           $mytext= $_POST[‘mytext'];
     $userid= $_POST['userid'];
     $newdir = “my_new_directory”;
     mkdir($newdir);
     $outf = $newdir."/mytext.txt";
     $falf=fopen($outf,"x");
     fprintf($falf,"%s",$mytext);
     fclose($falf);

You can also browse and upload files
 <form enctype="multipart/form-data" action=“test2.php" method="post"><pre>
    Your ID  <input type="text" name="userid"/>
    <input type="hidden" name="MAX_FILE_SIZE" value="500000" />
    Choose a first file to upload: <input name=“f1file" type="file" />
    Choose a second file to upload: <input name=“f2file" type="file" />
    <input type="submit" value="go" />
   </pre></form>
And handle with this PHP
        $userid= $_POST['userid'];
     $newdir = “mynewdir”;
     mkdir($newdir);
     $fileok = 0;
     $target_path = $newdir."/myfile1,txt";
     if(move_uploaded_file($_FILES[‘f1file']['tmp_name'], $target_path)) {
      echo “File ".  basename( $_FILES[‘f1file']['name'])." uploaded<br>";
     } else{
       echo "There was an error uploading the first File, please try again!<br>";
       $fileok = 1;
     }
     $target_path = $newdir."/myfile2.pdb";
     if(move_uploaded_file($_FILES[‘f2file']['tmp_name'], $target_path)) {
      echo “File ".  basename( $_FILES[‘f2file']['name'])." uploaded<br>";
     } else{
       echo "There was an error uploading the second file, please try again!<br>";
       $fileok = 1;
     }
*/

## Integer
$int = 1;

## Floats
$float = 2.2;

## Arrays
### Multidimensional Arrays (like matrices)

### Associative Arrays (like dictionary in Python)


# Calculation
echo "\n------------------ Calculation ------------------\n";
## Basic Four Calculation
echo 1+2 . "\n";
echo 2.1-1 . "\n";
echo 2.3*3.2 . "\n";
echo 3.2/8 . "\n";
echo 2**3 . "\n";


# Conditional Structure
echo "\n------------------ Conditional Structure ------------------\n";
## If
if (1 > 0) {
    echo "condition1 is true" . "\n";
}
elseif (1 < 0) {
    echo "condition2 is true" . "\n";
}
else {
    echo "conditions above are all false" . "\n";
}

## Switch
$colour = "red";
switch ($colour) {
    case "red":
        echo "Your favorite color is red!";
        break;
    case "blue":
        echo "Your favorite color is blue!";
        break;
    case "green":
        echo "Your favorite color is green!";
        break;
    default:
        echo "Your favorite color is neither red, blue, nor green!";
}

## Logical Operators
### AND
if ((true && 1 > 0) and true) {
    echo "AND operator: `or` `&&`" . "\n";
}
### OR
if (true || false or true) {
    echo "OR operator: `or` `||`" . "\n";
}
### XOR
if (true xor 1 > 0) {
    echo "Double true is not for xor?" . "\n";
}
elseif (true xor false){
    echo "XOR need only one true" . "\n";
}


# Loops
echo "\n------------------ Loop ------------------\n";
$i = 0;

## While Loop
while ($i <= 10) {
    $i = $i + 1;
}

## Do while Loop
### Tips: Remember the ; after the while condition
$i = 10;
do {
    $i = $i - 1;
} while ($i >= 1);

## For Loop
### Tips: The for loop is specific for counter in PHP and while is more about condition.
$i = 2;
for ($counter = 1; $counter <= 5; $counter = $counter + 1) {
    $i = $i ** 2;
    echo $i . "\n";
}

## Foreach Loop
$colors = array("red", "green", "blue", "yellow");

foreach ($colors as $value) {
    echo "$value";
}

## Some conditional structure
### Tips: There is no `exit` method to exit the loop but you can specify the break layer to achieve similar effect.
$i = 1;
while (true) {
    if ($i <= 5) {
        echo $i . " <= 5 \n";
        $i = $i + 1;
        continue;
    }
    elseif ($i > 5 and $i <= 10) {
        echo $i . " > 5 \n";
        $i = $i ** 2;
    }
    else {
        echo $i . " > 10";
        break 1;
    }
}


# Functions
echo "\n------------------ Functions ------------------\n";
$var0="Time: ";
## Define functions
function foo ($var1) {
    return date("l F jS Y", $var1);
}
## Call functions
echo $var0.foo(time());


# Classes and Objects
echo "\n------------------ Classes and Objects ------------------\n";
## Define a class
class Fruit {
    // Properties
    public $name;
    public $colour;

    // Methods for the class
    function set_name($name) {
        $this -> name = $name;
    }
    function get_name() {
        return $this -> name;
    }
}






?>
</pre>
</body>
</html>