#!/bin/bash

# Parses Subversion repository log messages and outputs authors and their 
# respective number of commits. 
#
# Output can be sorted:
#
#	-a		sort output by author name
#	-c		sort output by number of commits
#
# Output can also be limited to commits from the last year:
#
#	-y		display only commits from the last 365 days
#
# Will run on first 'path' parameter, or on current working directory if no 
# path is supplied.
#
# TODO:
# - Add option to bin authors and commits by year
# - Convert to using XML output from SVN log (i.e. `svn log --xml`)

## Functions
#

function display_help {
	echo "usage: ${script_name} [option] [path]" >&2
	case "${uname_out}" in
		Linux*)
			display_help_linux
			;;
		Darwin*)
			display_help_mac
			;;
		*)
			;;
	esac
}

function display_help_mac {
	echo "	-h	display help"
	echo "	-a	sort output by author name"
	echo "	-c	sort output by number of commits"
	echo "	-y	display only commits from the last 365 days"
}

function display_help_linux {
	echo "	-h | --help		display help"
	echo "	-a | --authors	sort output by author name"
	echo "	-c | --commits	sort output by number of commits"
	echo "	-y | --year		display only commits from the last 365 days"
}

function print_authors {
	i=0
	until [[ $i -eq "${#authors[@]}" ]]; do
		printf "%-15s %i\n" "${authors[$i]}" "${commits[$i]}"
		(( i++ ))
	done
}

function print_authors_sorted_by_author {
	local i=0
	until [[ $i -eq "${#authors[@]}" ]]; do
		printf "%-15s %i\n" "${authors[$i]}" "${commits[$i]}"
		(( i++ ))
	done |
	sort -n -k1
}

function print_authors_sorted_by_commits {
	local i=0
	until [[ $i -eq "${#authors[@]}" ]]; do
		printf "%-15s %i\n" "${authors[$i]}" "${commits[$i]}"
		(( i++ ))
	done |
	sort -rn -k2
}

## Initialize conditional variables
#
one_year_ago=""
sort_by_author=0
sort_by_commits=0
today=""

# Retrieve script name
script_name=$(basename "$0")

## Retrieve OS
#
uname_out=$(uname -s)
case "${uname_out}" in
	Linux*)
		options=$(getopt -n ${script_name} -o hacy --long help,author,commits,year -- "$@")
		;;
	Darwin*)
		options=$(getopt hacy $*)
		;;
	*)
		printf "operating system %s not currently supported\n" "${uname_out}" >&2
		exit 1
		;;
esac

## Handle options
#
valid_options=$?
if [[ $valid_options -ne 0 ]]; then
	display_help
	exit 1
fi
#echo $options # Uncomment for debugging
set -- $options

while :
do
	case "$1" in
		-h | --help)
			display_help
			exit 0
			;;
		-a | --author)
			sort_by_author=1
			shift
			;;
		-c | --commits)
			sort_by_commits=1
			shift
			;;
		-y | --year)
			today=$(date +%Y-%m-%d)
			one_year_ago=$(date -v -1y -v +1d +%Y-%m-%d) # Date adjusted 365 days
			shift
			;;
		--)
			shift
			break
			;;
		*)
			break
			;;
	esac
done

## Retrieve path
#
# If no path was given, use current working directory; if one or more were 
# given, just use the first one.
#
path="."
num_args=${#@}
if [[ $num_args -gt 0 ]]; then
	path=$1
fi

authors=()
commits=()

# Get entire log once
if [[ $one_year_ago != "" ]]; then
	log=$(svn log -r "{${today}}:{${one_year_ago}}" "${path}")
else
	log=$(svn log "${path}")
fi

data_lines=$(echo "${log}" | egrep "r[0-9]+") # Get only lines from log that contain info we care about

while IFS= read line; do
	author=$(echo "${line}" | awk -F "|" '{print $2}' | tr -d ' ')

	# Uncomment if we want to treat blank author as "unknown"
	if [[ "$author" == "" ]]; then
		author="unknown"
	fi

	have_seen_author=0
	for (( i=0; i<${#authors[@]}; ++i )); do
		if [[ "${authors[$i]}" == "$author" ]]; then
			have_seen_author=1
			break
		fi
	done

	if [[ $have_seen_author -eq 0 ]]; then
		authors+=("${author}")
		commits+=(1)
	else
		(( commits[$i]++ ))
	fi
done <<< "${data_lines}"

## Print results (sorted, if requested)
#
# NOTE: If multiple sort options are specified, the order or precedence is 
#		specified as follows.
#
if [ $sort_by_commits -eq 1 ]; then
	print_authors_sorted_by_commits
elif [ $sort_by_author -eq 1 ]; then
	print_authors_sorted_by_author
else
	print_authors
fi
