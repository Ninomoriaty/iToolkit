Guideline Level

<table><tr><td  style="background: #C8D066 !important;"><center><font style="font-size:45px; color:black !important;"><b>
    FAQ - Strange Questions and Solutions
    </b></font></center></td></tr></table>

# Installation Problem

1. Internet question when using `install.packages()`

   ```r
   Warning: unable to access index for repository https://r-forge.r-project.org/src/contrib:
     internet routines cannot be loaded
   ```

   Solution:

   1. check the internet setting

      1. VPN or VPS should be checked corresponding firewall and port setting

   2. check the CRAN repository

      1. make sure that the mirror is usable and correct
      2. reset the CRAN repository by `options` or `Tools > Global Options > packages > primary repository `
      3. HTTP secure problem: uncheck the `Tools > Global Options > packages > Use sercure download method for HTTP` option.

   3. Conda env problem(Vertified)

      1. (Re)Install rstudio in conda environment (usr/lib is not the same as the env/lib)

         ```bash
         conda activate [envname]
         conda install -c r rstudio # it would provide the dependencies for the rstudio in this environment,including the internet probelm
         
         ```
      
      Tips: To specify the Internet problem, the lack of `curl` may be the real cause fo this problem 
      
      conda install curl
      
         ```
      
      2. 
         ```

2. `curl` defection

   ```bash
   # curl loss
   sudo apt install curl
   # curl dependnecy: libcurl4-openssl-dev and libssl-dev
   sudo apt install libcurl4-openssl-dev
   sudo apt install libssl-dev
   
   # check further error message of the command
   ```

3. 



# How to add library path

Lk: https://www.accelebrate.com/library/how-to-articles/r-rstudio-library

## Temporarily

```R
myPaths <- .libPaths()   # get the paths
myPaths <- c(myPaths[2], myPaths[1])  # switch them
.libPaths(myPaths)  # reassign them
```

## Permanently

Tips: Add this code section to the Rprofile or the Rprofile.site in the R bin library for your usage

```R
myPaths <- .libPaths()
myPaths2 <- c("/home/ninomoriaty/anaconda3/envs/SDAenv/lib/R/library", myPaths)
.libPaths(myPaths2)
```

