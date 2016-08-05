#!/bin/bash
sudo mkdir /store/census/pgsql
sudo chown postgres:postgres /store/census/pgsql
sudo chmod 775 /store/census/pgsql
#sudo mkdir /store/census/pgsql/data


sudo -u postgres /usr/lib/postgresql/9.5/bin/initdb --auth=trust --pgdata=/store/census/pgsql/data --encoding=UTF8

sudo -u postgres /usr/lib/postgresql/9.5/bin/pg_ctl -w -D /store/census/pgsql/data -o "-p 5433" start

sudo -u postgres psql --host=localhost --port=5433 --command="CREATE ROLE gcensus WITH LOGIN CREATEDB PASSWORD 'stuff'"

sudo -u postgres createdb --host=localhost --port=5433 --owner=gcensus --encoding=UTF8 gcensus

sudo -u postgres psql --host=localhost --port=5433 --dbname=gcensus --command='CREATE EXTENSION postgis'

sudo -u postgres psql --port=5433 --dbname=gcensus --command='CREATE EXTENSION postgis_topology'

#sudo -u postgres psql --port=5433 --file=/usr/share/postgresql/9.5/contrib/postgis-2.2/topology.sql gcensus

# GRANT ALL PRIVILEGES ON SCHEMA topology to gcensus;
#  GRANT ALL PRIVILEGES ON sequence topology_id_seq to gcensus;
