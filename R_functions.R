###  All R functions:

## function to make individual groups
make_a_group <- function(neighb_list, gr_size, seed = 1L) {
  
  # if the seed is a 1 in character, make it integer to work after the first iteration
  if(seed=="1") seed <- 1L
  
  # make first group
  gr <- list(gr = c(as.numeric(names(neighb_list)[seed]), neighb_list[[seed]]),
             eval = as.numeric(names(neighb_list)[seed]))
  nn <- length(gr$gr)
  
  # Until you reach target cluster size, do:
  while(nn < gr_size){
    
    # find neighbours to evaluate
    to_eval <- gr$gr[!gr$gr %in% gr$eval]
    if(length(to_eval)==0) break() # if none, stop
    i = 1
    nb_to_eval <- length(to_eval)
    
    # and evaluate them sequentially until you did them al or reach target size
    while(nn < gr_size & i<=nb_to_eval){
      gr$gr <- unique(c(gr$gr, neighb_list[[as.character(to_eval[i])]]))
      i = i+1
      nn <- length(gr$gr)
    }
    
    # update the info on which polygon was evaluated
    gr$eval <- c(gr$eval, to_eval[1:(i-1)])
  }
  lapply(gr, sort)
}

## function that group the leftover polygons
##  (sometime, a polygon looses all it's neighbors because they were assign other groups)
group_empty_pol <- function(pol_name, neighb_list, the_groups){
  
  # find the original neighbors of this polygon
  all_ne <-neighb_list[[pol_name]]
  
  # Find the group which has the most polygons in common with its neighbors
  better_gr <- which.max(sapply(the_groups, function(xx) sum(xx %in% all_ne)))
  
  # add him to the group
  the_groups[[better_gr]] <- sort(c(as.numeric(pol_name), the_groups[[better_gr]]))
  the_groups
}


## Function to split the whole shape in groups from neighbours list

cluster_neighbours <- function(neighb_list, gr_size, seed = 1L){
  
  # validating the value of gr_size:
  if(gr_size<=1) stop("gr_size should be higher than 1")
  if(gr_size>=length(neighb_list)) stop("gr_size should be lower than the maximum number of polygons")
  
  
  # adding a seed to the code, 
  # the seed is the first polygon to be grouped.  It's an index for the table
  # as to be a integer from 1:n or "random" 
  if(length(seed)!=1) stop("seed should be a vector of size 1")
  if(class(seed) == "character") {
    if(seed!="random") 
      stop("The seed must be a integer (index of the wanted polygon) from 1:n or 'random'") 
  } else {
    if(!(1L<=as.integer(seed) & as.integer(seed)<=length(neighb_list)))
      stop("The seed must be a integer (index of the wanted polygon) from 1:n or 'random'")
  }
  
  # initialize object
  ll_group <- list()
  neighb_temp <- neighb_list
  
  # make island with only one polygon each their own group
  if(sum(sapply(neighb_temp, length)==0)>0) {
    ll_group <- neighb_temp[sapply(neighb_temp, length)==0] %>%
      names %>% 
      as.numeric() %>% 
      as.list
    neighb_temp <- neighb_temp[sapply(neighb_temp, length)!=0]
  }
  
  # validate that the seed chosen wasn't an island
  if(seed!=1 & seed!="random") if(!any(names(neighb_temp)==as.character(seed))) stop("The seed chosen is a island with no neighbors")
  
  # chose seed if random
  if(seed=="random") seed <- sample(names(neighb_temp), 1)
  
  # making the seed character to be extracted by name
  seed <- as.character(as.integer(seed))
  
  # until all are neighbors list is empty, do:
  while(length(neighb_temp)>0){
    
    # Make a group
    res1 <- make_a_group(neighb_list = neighb_temp, gr_size = gr_size, seed = seed)
    
    # reset seed to 1L so subsequent grouping follow previous group
    seed <- 1L
    
    # save the group
    ll_group <- append(ll_group, list(res1$gr))
    
    # Remove from neighbor list all polygons from the just made group
    neighb_temp <- neighb_temp[!names(neighb_temp)%in%as.character(res1$gr)] %>% 
      lapply(function(xx) xx[!xx%in%res1$gr])                                       # remove the used neighbors from neighbor list of left neighbors
    
    ## managing polygons with no more neighbors
    lost_all_neighb <- neighb_temp[sapply(neighb_temp, length)==0] %>% names
    
    # if some polygons lost all their neighbors, do:
    if(length(lost_all_neighb)>0){
      for(i in 1:length(lost_all_neighb)){
        #assign the lonely polygon to a group
        ll_group <- group_empty_pol(pol_name = lost_all_neighb[i], neighb_list = neighb_list, the_groups = ll_group)
        
        # Remove it from the neighbor list
        neighb_temp <- neighb_temp[!names(neighb_temp)%in%lost_all_neighb[i]] %>% 
          lapply(function(xx) xx[!xx%in%as.numeric(lost_all_neighb[i])])       
      }
    }
    ##
  }
  ll_group
}


## final grouping of tiny groups
##  (sometime, some tiny groups are left, we want to fuse them to bigger ones)
fuse_tiny_group <- function(initial_gr, init_sh, hard_min_size){
  
  # initialize data
  initial_gr <- initial_gr %>% 
    setNames(1:length(.))
  
  # assign good group to polygon, and dissolve the shape by group
  sh_with_gr <- bind_cols(init_sh,
                          do.call("rbind", mapply(function(gr, index) cbind(index, gr),1:length(initial_gr), initial_gr)) %>% 
                            as.data.frame() %>% 
                            arrange(index)
  ) %>% 
    select(gr) %>% 
    mutate(gr = as.character(gr)) %>% 
    group_by(gr) %>% 
    summarise(n = n()) %>% 
    arrange(as.numeric(gr)) 
  
  # Calculate new neighbors
  adj_gr <- sh_with_gr %>% 
    st_touches() %>% 
    setNames(sh_with_gr$gr)
  
  # grouping the island to the closest group distance wise
  dist_mat <- st_distance(sh_with_gr[sapply(adj_gr, length)==0,], sh_with_gr[sapply(adj_gr, length)!=0,])
  island_paired <- apply(dist_mat, 1, which.min) %>%
    magrittr::extract(sh_with_gr[sapply(adj_gr, length)!=0,]$gr, .) %>% 
    setNames(sh_with_gr[sapply(adj_gr, length)==0,]$gr)
  
  for(i in names(island_paired)){
    initial_gr[[island_paired[i]]] <- sort(c(initial_gr[[island_paired[i]]],initial_gr[[i]]))
    initial_gr[[i]] <- NULL
  }
  
  # find the small groups that are not island (which we just grouped)
  small_gr <- sh_with_gr$gr[sh_with_gr$n<hard_min_size]
  not_island <- names(adj_gr)[sapply(adj_gr, length)!=0]
  small_gr <- intersect(small_gr, not_island)
  
  # for every small group, fuse it with the smaller number of polygons
  for(igr in small_gr) {
    groups_close_by <- adj_gr[[igr]]
    
    group_selected <- initial_gr %>% 
      subset(.,names(.)%in%groups_close_by) %>% 
      sapply(length) %>% 
      sort() %>%
      names %>% 
      head(1)
    
    initial_gr[[group_selected]] <- sort(c(initial_gr[[group_selected]], initial_gr[[igr]]))
    initial_gr[[igr]] <- NULL
  }
  initial_gr
}
