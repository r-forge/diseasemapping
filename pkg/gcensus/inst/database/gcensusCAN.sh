#!/bin/bash


RPATH="/store/census/CAN"
for year in $RPATH/[[:digit:]]*; do
  shortYear="${year/$RPATH\//}";
  psql --username=gcensus --dbname=gcensus --host=localhost --port=5433 --command=CREATE\ SCHEMA\ can$shortYear
  for level in $year/[[:digit:]]*; do
    levelShort="${level/$RPATH\//can}";
    levelShort="${levelShort/\//.l}";
    echo $levelShort;
    #add boundary files
    for sfile in $level/*.shp; do
		psql --host=localhost --port=5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"map"
      shp2pgsql -s 4326 -W "LATIN1" $sfile $levelShort"map" gcensus | psql --quiet --host=localhost --port=5433 --username=gcensus --dbname=gcensus
    done

    #add ecumene files
    for sfile in $level/ecumene/*.shp; do
		psql --host=localhost --port=5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"ecumene"
        shp2pgsql -s 4326 -W "LATIN1" $sfile $levelShort"ecumene" gcensus | psql --quiet --host=localhost --port=5433 --username=gcensus --dbname=gcensus
    done

    #add population and income data if exists
    if [ -e "$level/population.dbf" ]
    then
		psql --host=localhost --port=5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"pop"
        shp2pgsql -s 4326 -W "LATIN1" -n $level/population.dbf $levelShort"pop" gcensus | psql --quiet --host=localhost --port=5433 --username=gcensus --dbname=gcensus
    fi

    if [ -e "$level/income.dbf" ]
    then
		psql --host=localhost --port=5433 --username=gcensus --dbname=gcensus --command=DROP\ TABLE\ $levelShort"inc"
        shp2pgsql -s 4326 -W "LATIN1" -n $level/income.dbf $levelShort"inc" gcensus | psql --quiet --host=localhost --port=5433 --username=gcensus --dbname=gcensus
    fi

  done
done

