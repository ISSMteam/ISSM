#/bin/bash
#Launch cron as a windows service, using cygrunsrv: 

cygrunsrv --install cron --path /usr/sbin/cron --args -n
net start cron
