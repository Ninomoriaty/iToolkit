#!/bin/bash

# su
## F: assume the identity of another user and either start a new shell session with that user's ID, or to issue a single command as that user. 
#$ su [options] [-] [user_name] 
# and then type the passwd would login in the administrator account
## Options
#-l, --login # login shell and change the home dir to the user in the argument. Default Superuser.
#-c, --command '<command>' # exec the command with the permissions of specified user. Tips: enclose the command in quotes

# sudo
## F: set up a configuration file called /etc/sudoers and define specific commands that particular users are permitted to execute as a different user (usually the superuser) under an assumed identity.
## Ab: 
## 1. sudo does not require access to the superuser's password, but requires the user’s own password.
## 2. does not start a new shell, nor does it load another user's environment
#$ sudo [commands]  # use the administrator authority for particular commands
## Options
#-i  # start an interactive superuser session
#-l  # see the priviledge of sudo
