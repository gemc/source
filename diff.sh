#!/bin/zsh

# script to show differences between clas21Tags/source and gemc/source
#
# if the prompt argument is given, then the script will ask for user input
# to copy each file into the tag
#
prompt="no"

if [[ $# -gt 0 ]]; then
	if [[ $1 == "prompt" ]]; then
		prompt="yes"
	fi
fi

ignores="-x .idea -x .git -x .gitignore -x *.o -x moc_*.cc -x *.a -x api -x .sconsign.dblite -x diff.sh "

printf "Ignoring $yellow$ignores$reset\n\n"

diffs=$(diff -rq $=ignores ../clas12Tags/source . | sed 's/Files //g' | sed 's/ and / /g' |  sed 's/ differ//g')

# create an array from diffs where the discriminator is carriage return
diffs=("${(@f)diffs}")

print "\nDiffs:\n"
for d in $diffs; do
	source=$(echo "$d" | awk '{print $1}')
	target=$(echo "$d" | awk '{print $2}')
	if [[ $prompt == "yes" ]]; then
		clear
		printf "\nDiffs of source: $yellow$source$reset with $yellow$target$reset:\n"
		diff $source $target
		printf "\n$magenta Copy? (y/n)$reset\n"
		read -r answer
		echo $answer
		if [[ $answer == "y" ]]; then
			cp $source $target
		fi
	else
		printf "$d\n"
	fi
done


printf "\n- Changing initializeBMTConstants, initializeFMTConstants, initializeRTPCConstants to initialize before processID"
sed -i s/'initializeBMTConstants(1)'/'initializeBMTConstants(-1)'/ hitprocess/clas12/micromegas/BMT_hitprocess.cc
sed -i s/'initializeFMTConstants(1)'/'initializeFMTConstants(-1)'/ hitprocess/clas12/micromegas/FMT_hitprocess.cc
sed -i s/'initializeRTPCConstants(11)'/'initializeRTPCConstants(-1)'/ hitprocess/clas12/rtpc_hitprocess.cc
