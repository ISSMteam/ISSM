#!/bin/sh
#
# SOURCE_SEEK.SH - Look for IP addresses and mention of 'password' in source files - 2009-10-05
#
if [ -n "${1}" ]; then
  clear
  echo "Checking for IP addresses"
  echo "--------------------------------------------"
  find $1 -name "*" -exec grep -il "128\.149\." {} \;
  find $1 -name "*" -exec grep -il "137\.78\."  {} \;
  find $1 -name "*" -exec grep -il "137\.79\."  {} \;
  find $1 -name "*" -exec grep -il "137\.228\." {} \;
  echo "--------------------------------------------"
  echo "\nChecking for 'password'"
  echo "--------------------------------------------"
  find $1 -name "*" -exec grep -il "password"   {} \;
  echo "--------------------------------------------"
  echo "\nChecking for JPL user IDs"
  echo "--------------------------------------------"
  find $1 -name "*" -exec grep -ilf jplusers.db {} \;
  echo "--------------------------------------------"
else
  echo -e "\nUSE: $0 <directory_of_source_files>\n";
fi
#
# EOF - 2011-08-05 - Jay Dyson <jdyson@jpl.nasa.gov>
