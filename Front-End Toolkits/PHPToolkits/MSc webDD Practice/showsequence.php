<?php

echo <<<_HEAD

<head>
   <meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Sequence Display</title>
<link href="style.css" rel="stylesheet" type="text/css" />
</head>

<body>
<img border="0" hspace="0" src="titlebar.gif" width="800" height="200">

_HEAD;

if(isset($_POST['pdbid']) && isset($_POST['chnid']))
{
    $pdbid = $_POST['pdbid'];
    $chnid = $_POST['chnid'];
    printf("<pre>pdbid is %s and chainid is %s\n",$pdbid,$chnid);
    $lpdb = strlen($pdbid);
    $lchn = strlen($chnid);
    if(($lpdb != 4) && ($lchn != 1))
    {
        printf("pdbid must be 4 characters(now %d), chain id must be one character(now %d)",$lpdb,$lchn);
    } else {
        $pathel = substr($pdbid,1,2);
        $lpathel = strtolower($pathel);
        $lpdbid = strtolower($pdbid);
        $comtodo = "zcat pdb/$lpathel/pdb${lpdbid}.ent.gz | /usr/progs/pdbcheck/pdbcheck_a $chnid";
        system($comtodo);
    }
    printf("</pre>");
}

echo <<<_EOP

   <form action="showsequence.php" method="post">
   <pre>
       PDB code<input type="text" name="pdbid"/>
       Chain Id<input type="text" name="chnid"/>
                   <input type="submit" value="go" />
   </pre>
   </form>

_EOP;

echo <<<_TAIL
</body>
</html>
_TAIL;

?>



