# Thermal-storage-tank
Thermal storage tank simulating transient behavior including thermal stratification.















## Tools Used

 - [TMUX- terminal multiplexer](https://hackernoon.com/a-gentle-introduction-to-tmux-8d784c404340)
 - [VirtualEnv- virtual environment](https://virtualenv.pypa.io/en/latest/userguide/) 


## Basic Linux Commands
- `$ ls` list all files in the current directory
- `$ cd subfolder` changes your current directory to subfolder
- `$ cd ../` changes your current directory to the one above it
- `$ cd ~` changes your current directory to be the home directory
- `$ mkdir folder` create a folder in the current directory named folder
- `$ rm file.txt` deletes file.txt
- `$ rm -rf directory` deletes directory and anything inside it
- `$ mv src tgt/` moves directory src to tgt directory. 
- `$ mv old new` renames directory old to new. 
- `$ chmod 777 script.sh` changes permissions for script.sh to be Read, Write, and Executable for all users
- `ctrl+c` cancels the currently running command (stops a program)
- `$ exit` exits and closes the current ssh session


## pip3
` $ sudo apt install python3-pip`
pip is the python package manager that we use. We must be careful about which version of python pip is installing for. Be sure to check the pip version if you are unsure.
- `$ pip install package` Installs a package from PyPI
- `$ pip uninstall package` Uninstalls a package
- `$ pip list` Lists all install packages
- `$ pip --version` Prints out the version
- `$ pip install -r requirements.txt` Installs all packages listed in requirements.txt
- `$ pip freeze` List out all installed packages and their version numbers


## tmux

Tmux is a terminal multiplexer, meaning that it allow us to have multiple terminal windows running at the same time. This is useful for running programs that output to the command line, or for running programs continuously. This way, when the ssh session is closed, the process does not die.

The general workflow to start such a program is to create a new tmux session, run the program from that session, and then detach from the session.

`$ tmux new -s session` create new session with name "session"

`$ tmux a -t session` attach to existing session with name "session"

`$ tmux ls` list all running sessions

`$ tmux kill-session -t session` deletes the session "session"

Once inside a tmux session
- `ctrl+b d` Detach session (return to main terminal)
- `ctrl+b [` Enter copy mode (to scroll up)
- `ctrl+b q` Quit/Exit copy mode

[tmux cheat sheet](https://tmuxcheatsheet.com/)

## virtualenv
A virtual environment provides an isolated development environment from which we can run code. When working inside a virtual environment, Python packages are only installed for that environment. This prevents package version and configuration conflicts.

`pip install virtualenv`

Creating and entering a virtualenv:

1. `$ virtualenv env` Creates a virtual environment named "env". This creates a folder with that name in your current directory, where all the python libraries and executables are stored for that environment.
2. `$ source env/bin/activate` Will enter the env virtual environment. If you have done it successfully, you will see the name of the environment at the beginning of your cursor line. 
3. Once you are inside a virtual environment, you can install or run programs as you do normally. 
4. `(env) $ deactivate` will exit the virtual environment, if you are in one.

Virtual environment only affects the commands and the packages that you install while inside it. Changes to files or directories will be kept after you exit the environment. It is good practice to do your work inside a virtual environment as that will keep your main Python install clean and less likely to have convoluted path errors.

## Vim
Vim is a text editor that runs from the command line. This means that we can use it to look at files and make modifications.

`$vim file.txt` opens file.txt in vim. The file does not exist, it will be created with that name.
Upon opening the file, you will be put into read-only mode. Use the arrow keys to move the cursor around.

To begin command mode, start with `esc`, then you can enter a command
- `i` enters insert mode. You can now freely write or edit the contexts of the text.
- `:q` exits the editor. q is short for quit
- `:q!` force exits the editor. Any changes made will not be saved.
- `:x` saves and exits the editor
- `:$` or `G` jumps to the end of the file
- `gg` or `:1` jumps to the start of the file

## Git
Git is a version management tool geared towards software development. Separate code branches can be worked on independently and then merged when they are completed. This means that developers can freely made modifications to a project without it affecting the live version. Git generally comes preinstalled on linux machines, but will need to be installed on [Windows](https://git-scm.com/download/win). Github is the most popular hosting service for developers to store their git repositories on.

### Initializing a repository
- `git clone http://github.com/user/repo.git` clones an existing repository into your local directory. Changes you make on that repository will be tracked by git
- [How to add a project to Github](https://help.github.com/en/articles/adding-an-existing-project-to-github-using-the-command-line)

The `.git` hidden folder in a directory signifies that the repository is being tracked by git. To stop a directory from being tracked, you can simply delete the .git folder. Other directories will be unaffected.

### Working in a repository
When working, it is best practice to make regular commits to your working branch.

Use `git status` to see the current status of the repository. It will tell you the current branch you are on, as well as the number of your commits in relation to the master branch.

When you are ready to make a checkpoint:
1. `git add .` adds all files in the directory and subdirectories to be tracked by git.
2. `git commit -m "message"` commits the current changes to the local branch, with a message. This is analogous to making a check point that can easily be reverted later to.
3. `git push` pushes the local branch to the associated one on the remote server. This will upload all changes made since the code was last branched, and attempt to merge the local and remote branch. This will ensure that the code is backed up on the remote server. Often times, many commits will be made before they are all pushed.

### Working on a branch
-  `git branch feature` creates a new branch named feature
- `git branch -a` prints out all the available branches
- `git checkout feature` switches to the feature branch. **Note the repository will change** to where the branch was when you last left it. Changes you make now, will be reflected in the feature branch.
