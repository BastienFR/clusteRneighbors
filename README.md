# clusteRneighbors
R functions to spatially cluster polygons that are touching each other aiming at a target cluster size.


# Example

```
library(sf)
library(dplyr)


## read the shapefile (mine from: http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lda_000b16a_e.zip)
DA <- st_read("C:/Users/DXD9163/Desktop/DA_regrouping/DA_CAN_2016_reduced.shp", stringsAsFactors = F, options = "ENCODING=windows-1252") %>% 
  filter(!is.na(DAUID))

## select only a part of it
DA_qc <- DA[grepl("Quebec",DA$PRNAME),]

## change from singlepart to multipart (polygon to multipolygon)
DA_qc <- DA_qc %>% 
  mutate(area = st_area(DA_qc)) %>% 
  group_by(DAUID) %>% 
  summarise(area=sum(area))

## create a neighbor list
neighb <- st_touches(DA_qc, sparse = T) %>% 
  setNames(1:nrow(.)) 

## make the first grouping with the required target size (it will be a minimum but for some exceptions see below))
gr1 <- cluster_neighbours(neighb_list = neighb, gr_size = 1000)

## evaluate the size of the groups
sapply(gr1, length)
## multiple groups are very small, these can be island or leftover polygons, we need to regroup them to larger group

## Fuse the tiny groups to bigger groups, you can specify a new minimum smaller than the first one, 
##  higher this minimum is, bigger can be your final groups (potentially ~ gr_size + hard_min_size)
final_group <- fuse_tiny_group(initial_gr = gr1, init_sh = DA_qc, hard_min_size = 400)

## the new sizes
sapply(final_group, length)


## Assign the proper group to every polygon and save the shape 
bind_cols(DA_qc,
          do.call("rbind", mapply(function(gr, index) cbind(index, gr),1:length(final_group), final_group)) %>% 
            as.data.frame() %>% 
            arrange(index)
) %>% 
  select(gr) %>% 
  mutate(gr = as.character(gr)) %>% 
  st_write("C:/Users/DXD9163/Desktop/DA_regrouping/sh_separe/result2.shp")
```


