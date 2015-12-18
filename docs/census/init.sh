#!/bin/bash
sudo mkdir /store/census/pgsql
sudo chown postgres:postgres /store/census/pgsql
sudo chmod 775 /store/census/pgsql
#sudo mkdir /store/census/pgsql/data


sudo -u postgres /usr/lib/postgresql/9.4/bin/initdb --auth=trust --pgdata=/store/census/pgsql/data --encoding=UTF8

sudo -u postgres /usr/lib/postgresql/9.4/bin/pg_ctl -w -D /store/census/pgsql/data -o "-p 5433" start

sudo -u postgres psql -p 5433 --command="CREATE ROLE gcensus WITH LOGIN CREATEDB PASSWORD 'stuff'"

sudo -u postgres createdb -p 5433 --owner=gcensus --encoding=UTF8 gcensus

sudo -u postgres psql -p 5433 --dbname=gcensus --command='CREATE EXTENSION postgis'

