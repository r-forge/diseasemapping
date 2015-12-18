#!/bin/bash




RPATH="/store/census/USA"
for year in $RPATH/[[:digit:]]*; do
  shortYear="${year/$RPATH\//}";
  psql --username=gcensus --dbname=gcensus -p 5433 --command=CREATE\ SCHEMA\ usa$shortYear
  for level in $year/[[:digit:]]*; do
	levelShort="${level/$RPATH\//USA}";
	levelShort="${levelShort/\//.l}";
  	for sfile in $level/*.shp; do
		echo $sfile;
		psql -p 5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"map"
	  shp2pgsql -I -s 4326 -W "LATIN1" $sfile $levelShort"map" gcensus | psql --quiet -p 5433 --username=gcensus --dbname=gcensus
	done
	if [ -e "$level/pop.dbf" ]
	then
		psql -p 5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"pop"
		shp2pgsql -n -W "LATIN1" $level/pop.dbf $levelShort"pop" gcensus | psql --quiet -p 5433 --username=gcensus --dbname=gcensus
	fi
	if [ -e "$level/income.dbf" ]
	then
		psql -p 5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"income"
		shp2pgsql -n -W "LATIN1" $level/income.dbf $levelShort"income" gcensus | psql --quiet -p 5433 --username=gcensus --dbname=gcensus
	fi
  done
done

