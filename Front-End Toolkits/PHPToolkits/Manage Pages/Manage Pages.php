<?php
# Session mechanism to share variables from different pages
session_start();
echo<<<_HEAD1
<html>
<body>
_HEAD1;

if(isset($_POST['fn']) && isset($_POST['sn'])) {
    # Inputs for sessions
    $_SESSION['forname'] = $_POST['fn'];
    $_SESSION['surname'] = $_POST['sn'];
}
elseif (!(isset($_SESSION['forname']) && isset($_SESSION['surname'])))
{
    # Reminder for mistaken inputs
    # header(); is the function to
    header('xxxx.php');
}

# Empty the sessions
$_SESSION = array();
if( session_id() != "" || isset($_COOKIE[session_name()]))
    setcookie(session_name(), '', time() - 2592000, '/');
session_destroy();


?>