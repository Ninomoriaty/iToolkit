<?php
$db_hostname = 'localhost';
$db_database = 'mytest';
$db_username = 'prt';
$db_password = 'tef7glu';
$db_server = mysql_connect($db_hostname,$db_username,$db_password);
if(!$db_server) die("Unable to connect to database: " . mysql_error());
mysql_select_db($db_database,$db_server) or die ("Unable to select database: " . mysql_error());
$query = "select * from students";
$result = mysql_query($query);
if(!$result) die("unable to process query: " . mysql_error());
$rows = mysql_num_rows($result);
for($j = 0 ; $j < $rows ; ++$j){
    $row = mysql_fetch_row($result);
    $sn = $j + 1;
    echo "Student $sn $row[2] $row[1] \n";
}

mysql_close($db_server);

