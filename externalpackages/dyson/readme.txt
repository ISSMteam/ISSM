-- ARCHIVE MANIFEST ------------------------------------------ 2011-08-19 --

  Filename        Description
  --------------  ----------------------------------------------------------
  readme.txt      This file: describes script functions.
  ldap.pl         PERL script for acquiring JPL user IDs via cron job.
  source_seek.sh  Code-scanning script.  Surveys scripts for hard-coded IP
                  addresses, password callouts, and user IDs.

-- SCRIPT INFORMATION ------------------------------------------------------

  LDAP.PL (requires PERL, UNIX ldapsearch, and related utilities)
  --------------------------------------------------------------------------
  The LDAP.PL script is intended solely for twice-daily invocation via cron.
  This script should be run once at noon, once at midnight.  This script 
  generates the jplusers.db file on which the SOURCE_SEEK.SH shell script 
  relies for its scanning of source files.

  SOURCE_SEEK.SH (requires UNIX shell and related utilities)
  --------------------------------------------------------------------------
  The SOURCE_SEEK.SH shell script is intended for regular cron invocation,
  but may be run at system administrator discretion.  This script uses UNIX
  'find' and 'grep' utilities to scan code files for hard-coded JPL Internet
  Protocol (IP) addresses, password callouts, and JPL user IDs.  This script
  accepts the path where source code files reside as a command line
  argument (e.g., ./source_seek.sh /path/to/sourcecode/files).

-- AUTHOR CONTACT INFO -----------------------------------------------------

  Please direct questions regarding these files to the developer.
 
  Name   : Jay Dyson, CISSP
  Role   : IT Security Engineer
  Group  : JPL IT Security
  E-mail : jdyson@jpl.nasa.gov
  Phone  : 818-397-4960

------------------------------------------------ LAST UPDATED: 2011-08-19 --
