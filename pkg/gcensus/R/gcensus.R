#' Retrieve and combine us boundary, population and income data

#' @param country string specifying a country
#' @param year census year
#' @param level see below
#' @param subset passed as \code{where} to postgreql
#' @param spat spatial extent to crop to
#' @param database details of database to retrieve data from
#' @param income logical, retrieve income data
#' @details 

#' level	abbreviation	name
#' 0	USA	United States of America
#' 1	STATE	State
#' 1.1	PUMA	Public Use Microdata Areas
#' 1.2	ZCTA	Zip Code Tabulation Areas
#' 2	COUNTY	County
#' 2.1	COSUB	County Subdivision
#' 3	TRACT	Tract
#' 4	BLKGRP	Block Group


#' @export
gcensus = function(
		country='USA',
		year=2010, 
		level=1,
		subset = NULL, 
		spat = NULL,
		database = c(dbname='gcensus', user='gcensus', host='localhost',  port=5433),
		income=TRUE,
		t_srs=NULL,
		simplifyTolerance = 0) {
  
	level = gsub("\\.", "_", level)
	idString = paste('id', level, sep='')
	
	tempFileShp = paste(tempfile(), '.shp', sep='')
	
	spat_srs = NULL
	if(!is.null(spat)) {
		spat_srs = projection(spat)
		if(is.na(spat_srs))
			spat_srs = crsLL
		if(is.null(t_srs))
			t_srs = spat_srs
		spat = extent(spat)
		spat = c(xmin(spat), ymin(spat), xmax(spat), ymax(spat))
	}
	
	if(length(subset)) {
	  where = paste(names(subset), " like '", unlist(subset), "'", collapse="AND", sep='')
		whereFull = paste(
				'where', where
		)
	} else {
		where = whereFull = NULL
	}
	
	
	layer = paste(tolower(country), year, ".l", level, sep="")
	layerMap = paste(layer, "map",sep="")
	
	
  if(length(database)){
		
		forOgr = list(
				src_datasource_name=
 	    			paste("PG:", 
 	      				paste(names(database), '=', database, collapse=' ')),
 	  		dst_datasource_name=tempFileShp,
  			layer=layerMap,
				simplify=simplifyTolerance
  	)
		for(D in c('t_srs', 'where','spat','spat_srs')) {
			Dhere = get(D)
			if(!is.null(Dhere))
				forOgr[[D]] = Dhere
		}
		
		do.call(
				gdalUtils::ogr2ogr,
				forOgr
		)
		
 	  
		map=rgdal::readOGR(dirname(tempFileShp), 
    		gsub("\\.shp$", "", basename(tempFileShp)), 
    		stringsAsFactors=FALSE,verbose=FALSE)
		# overwrite projection because sometimes it's wrong
		if(!is.null(t_srs)) projection(map) = crs(t_srs)
  	
		if(is.null(whereFull)) {
			# should be 'id' instead of idString if the database is clean
			whereFull = paste("WHERE ", idString, " IN (", 
					paste(unique(map@data[['id']]), collapse=','), ")",
					sep='')
		}
		
  	con <- do.call(
				DBI::dbConnect,
				c(
						list(drv=RPostgreSQL::PostgreSQL()),
						as.list(database)
				)
		)
  	
		res <- DBI::dbSendQuery(con, 
    		paste("select * from ", tolower(country), year, ".l", level, 
    				"pop ", whereFull, 
						sep='')
		)
  	
  	pop = DBI::fetch(res, -1)
		
	  
    inc2 = pop[,grep("gid", names(pop), invert=TRUE)]
  	if(income) {
    	res <- try(DBI::dbSendQuery(con, 
      				paste("select * from ", tolower(country), year, 
        					".l", level, "income ", whereFull, sep='')
    			))
    	if(class(res)!= 'try-error') {
    		inc = DBI::fetch(res, -1)
    		inc2 = merge(inc[,-1], pop[,-1])
			}
  	}
		
  	DBI::dbDisconnect(con)
		
  } else {
	  warning("haven't written this yet")
  }
	
	if(idString %in% names(inc2))
		inc2$id = inc2[[idString]]
	
	inc2$id = as.character(inc2$id)
	map$id = as.character(map$id)
	
  map2=sp::merge(map, inc2, by='id')
	
	rownames(map2@data)= as.character(map2$id)
	map2 = sp::spChFIDs(map2, rownames(map2@data))
	
  map2
}	

if(FALSE) {	
	# simplify topology
	# see https://trac.osgeo.org/postgis/wiki/UsersWikiSimplifyWithTopologyExt
#	Steps ¶
	
#	Creates the target table that will contain simplified geometry. Initial departement polygons are dumped to handle simple objects.
	temptable = 'tempsimplifytopology'
	tempsimplified = 'foo'
	res <- try(DBI::dbSendQuery(con, 
    			paste('drop table ', temptable, ';',
							sep='') ), TRUE)
	
	res <- DBI::dbSendQuery(con, 
  		paste('create table ',
					temptable, ' as (select * from ', 
					layerMap, ' ', whereFull, '); ',
					'create index temp_geom_gist on ', temptable, ' using gist(geom);',
					sep='') )
	
#	-- adds the new geom column that will contain simplified geoms
	res <- DBI::dbSendQuery(con, 
			paste('alter table ', temptable, 
					" add column simple_geom geometry(POLYGON, '4326');",
					sep='')
	)
	
#	Creates a topology from the departements
# -- create new empty topology structure
	res <- try(DBI::dbSendQuery(con, 
					"select topology.DropTopology('topok');"
			))		
	stuff = try(DBI::fetch(res, -1))
	
	
	res <- DBI::dbSendQuery(con, 
			paste("select topology.CreateTopology('topok',4326,",
					simplifyTolerance, ");", sep=''
			))
	stuff = DBI::fetch(res, -1)
	
	
#	-- add all departements polygons to topology in one operation as a collection
	# this takes a long time
	res <- DBI::dbSendQuery(con, 	
			paste("select ST_CreateTopoGeo('topok',ST_Collect(geom))	from ", temptable,";", sep='')
	)
	stuff = DBI::fetch(res, -1)
	
#	Create a new topology based on the simplification of existing one. (should not be the right way to do it, but calling ST_ChangeEdgeGeom)
	res <- try(DBI::dbSendQuery(con, 
					"select topology.DropTopology('topok2');"
			))
	stuff = try(DBI::fetch(res, -1))
	
	res <- DBI::dbSendQuery(con, 	
			paste("select CreateTopology('topok2',4326,0);",
					"select ST_CreateTopoGeo('topok2', geom)",
					"from (",
					"select ST_Collect(st_simplifyPreserveTopology(geom, 10000)) as geom",
 					"from topok.edge_data) as  ", tempsimplified,";", sep=''
			)
	)
	stuff = DBI::fetch(res, -1)
	
# Retrieves polygons by comparing surfaces (pip is not enough for odd-shaped polygons)
	"with simple_face as (
  		select st_getFaceGeometry('topok2', face_id) as geom
			from topok2.face
			where face_id > 0
			) update new_dept d set simple_geom = sf.geom
			from simple_face sf
			where st_intersects(d.geom, sf.geom)
			and st_area(st_intersection(sf.geom, d.geom))/st_area(sf.geom) > 0.5;"
} # end simplify topology



