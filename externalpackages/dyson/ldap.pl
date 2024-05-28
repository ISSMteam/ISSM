#!/usr/bin/perl
#
# LDAP.PL - Acquire JPL LDAP data and pare down to user IDs - 2011-08-01
#
$cat="/bin/cat"; $srt="/bin/sort"; $ldb="/usr/bin/ldapsearch";                 # Binary locations (change as needed)
$lds="ldap.jpl.nasa.gov"; $ou="ou=personnel,dc=dir,dc=jpl,dc=nasa,dc=gov"; # Define LDAP server & organization
$bn="jplbadgenumber"; $txt="$bn".".txt"; $db="jplusers.db";                # Initialize values
#
system(`$ldb -x -h $lds -b $ou uid=* uid > $txt`);                         # Execute LDAPsearch, write to file
#
open(I,"<$txt");                                                           # Open input file handle
  open(O,">badges.tmp");                                                   # Open output file handle
    while (<I>) { chomp($_);                                               # Spool through input file
      if ($_ =~ "^uid: ") {                                                # If line begins with UID value,
        $_=~s/uid: //; $_=~s/,(.*)$//; $_=~s/\n//;                         #   strip line of extraneous data
        if (!($_=~/^\d\d\d\d\d\d/)) { print O "$_\n"; }                    # If UID is non-numeric, print output
      }                                                                    # Close condition
    }                                                                      # Close spool
  close(O);                                                                # Close output file handle
close(I);                                                                  # Close input file handle
#
system(`$cat badges.tmp |$srt -u > $db 2>/dev/null`);                      # Sort output
unlink("badges.tmp"); unlink("$txt");                                      # Unlink temp files
#
# EOF - 2011-08-01 - Jay Dyson <jdyson@jpl.nasa.gov>
