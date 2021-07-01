<?php
include "./start1.php";
if(isset($_POST['natmax'])) {
    echo "The value set for NATMAX in the form was ".$_POST['natmax'];
}
else {
    echo "There was no value for NATMAX set in the form";
}
?>