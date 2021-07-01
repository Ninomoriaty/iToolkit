<?php

# echo <<<_Name is a method to embed HTML within PHP files, which begins with the echo <<<_Name and ends with _Name
# Head
echo <<<_HEAD

<html>
<head>
# Actually there should be some metadata for HTML file
</head>
<body>

_HEAD;

# The main body of PHP
if(isset($_POST['yname'])) {
    # Response to the inputs
    echo "<p>You typed ",$_POST['yname']," in the box</p>";
    # Restart the main page and as the input is not accepted again, Go to the main page.
    echo '<a href="basic8.php">go back to the original page</a>';
}
else {
    # The Main Page wait for inputs
    echo <<<_FORM
    # Form needed inputs
    <form action="basic8.php" method="post">
    <pre>
           your name<input type="text" name="yname"/>
           <input type="submit" value="go" />
    </pre>
    </form>

    _FORM;
}

# Tail
echo <<<_TAIL

</body>
</html>

_TAIL;

?>