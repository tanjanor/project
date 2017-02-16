#create a file structure for a project

#scripts: directory for all python/R scripts
#bash: driver scripts that call all other scripts and execute pipelines
#input: input directory
#output: output directory
#logs: stdout and stderr of runall.sh scripts
#datasets: directory for storing all datasets
#general_scripts: directory for storing scripts that are use in many projects
#bin: all small scripts that should be in PATH


mkdir "$1"
cd "$1"
mkdir scripts bash input output logs
mkdir datasets general_scripts bin
> readme.txt
cd ./bash/
> runall.sh

