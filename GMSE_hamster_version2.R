library(GMSE); # Note that the semi-colons are not really necessary

#===============================================================================
# BD: Let's simplify by defining the different months explicitly in the code
TMAX <- 24; # Maximum number of months in the simulation;
JAN  <- seq(from = 1,  to = TMAX, by = 12); # Time steps in JAN, FEB, MAR, etc.
FEB  <- seq(from = 2,  to = TMAX, by = 12);
MAR  <- seq(from = 3,  to = TMAX, by = 12);
APR  <- seq(from = 4,  to = TMAX, by = 12);
MAY  <- seq(from = 5,  to = TMAX, by = 12);
JUN  <- seq(from = 6,  to = TMAX, by = 12);
JUL  <- seq(from = 7,  to = TMAX, by = 12);
AUG  <- seq(from = 8,  to = TMAX, by = 12);
SEP  <- seq(from = 9,  to = TMAX, by = 12);
OCT  <- seq(from = 10, to = TMAX, by = 12);
NOV  <- seq(from = 11, to = TMAX, by = 12);
DEC  <- seq(from = 12, to = TMAX, by = 12);
# We can now use each month in the loop to identify the month
#===============================================================================

#===============================================================================
# BD: GMSE doesn't do carrying capacity by cell, but we can force it to here
#     by creating a new function to deal with the resource_array
#     This should probably be only run after movement has happened and new nests
#     are built so that offspring are not removed before reproductively active
#     This function is quite slow -- probably want to speed it up somehow
cell_K <- function(res, DIM_1 = 447, DIM_2 = 447){
  for(i in 1:DIM_1){
    for(j in 1:DIM_2){
      sum_ij <- sum(res[, 5] == i & res[, 6] == j)
      if(sum_ij > 1){
        on_ij             <- which(res[, 5] == i & res[, 6] == j);
        loser_ij          <- sample(x = on_ij, size = sum_ij - 1);
        res[loser_ij, 7]  <- 0; # Stop in tracks
        res[loser_ij, 9]  <- 1; # Definite death
        res[loser_ij, 10] <- 0; # No birth
        res[loser_ij, 11] <- 0; # No birth
      }
    }
  }
  return(res);
}
#===============================================================================

#===============================================================================
# BD: Add crop type to the top landscape layer. The 'land' in GMSE is a 3D array
#     with the farm number in layer 3, the crop amount in layer 2, and nothing 
#     in layer 1. The code below redefines layer 1 by putting a crop type (1 or
#     2) down for each farm number (you can change these as desired).
#     There are better, more general, ways to code this, but basically the 
#     function below puts crop type 1 for farms 1 and 3, and crop type 2 for 
#     farms 2 and 4. You can add more farms and crops with more `if` statements.
#     See below that the function is run before looping over time.


#Imke: I changed the crop type function so it now assigns a random crop type to each farm at the start of the simulations

crop_type <- function(land){
  xdim <- dim(land)[1];
  ydim <- dim(land)[2];
  for(i in 1:xdim){   
  for(j in 1:ydim){
        if (land[i, j, 3] %in% 1:4) {  #1:number of farms/stakeholders 
          land[i, j, 1] <- sample(1:3, 1) #1:number of possible crop types 
      }
    }
  } 
    return(land);
} 
#===============================================================================

#===============================================================================
# BD: Have something happen to hamsters on a particular cell crop type. The code
#     below checks each hamster's cell location, then sees if the crop on that
#     location is a 1 or a 2 (feel free to add more `if` statements), then it
#     adjusts the birthrate (column 10) depending on the crop type. You can also
#     edit other things, of course, like death probability (column 9). I show 
#     how to use the function in the month of July below.

#Imke: changed parameters to simulate crop types 1-3
#I needed to make two functions, one to change parameters for types 1 and 2 in a monthly basis (e.g. no survival at all in type 1 fields)

crop_hamster_1 <- function(land, res){     
    N <- dim(res)[1]; # Number of hamsters
    for(i in 1:N) { 
        xloc <- res[i, 5];
        yloc <- res[i, 6];
        crop <- land[xloc, yloc, 1];
        
        cat("Hamster ", i, " at (", xloc, ", ", yloc, ") has crop type: ", crop, "\n")
        
          if (crop == 1) {
          res[i, 9] <- 1  # 100% mortality
        } else if (crop == 2) {
          res[i, 9] <- res[i, 9] + 0.10  # Increase mortality by 10%
        } else {
          #default case, nothing changes
      }
    }
    return(res);
}



#changes in birth rate only happen in certain months, this function is designed for those months: 

crop_hamster_2 <- function(land, res) {
  N <- dim(res)[1]  # Number of hamsters
  for (i in 1:N) {   
    xloc <- res[i, 5]
    yloc <- res[i, 6]
    crop <- land[xloc, yloc, 1]
    
    cat("Hamster ", i, " at (", xloc, ", ", yloc, ") has crop type: ", crop, "\n")
  
    if (crop == 1) {
      res[i, 9] <- 1  # 100% mortality 
    } else if (crop == 2) {
      res[i, 9] <- res[i, 9] + 0.10  # Increase mortality by 10%
    } else if (crop == 3) {
      res[i, 10] <- 1.18  # Set birthrate to 1.18 (extra litter)
    } else {
      #default
    }
  }
  return(res)
}

#===============================================================================
#parameters specified here are for crop type 3, except the extra litter which is modelled by the function crop_hamster_2
#for types 1 and 2, a different mortality is modelled each month (also look into no births in type 1 and different mortality in type 3)

DIM_1 <- 447; # Land dimension 1
DIM_2 <- 447; # Land dimension 2
# BD: Initialise the first output
sim_old   <- gmse_apply(stakeholders    = 4,
                        res_movement    = 0, 
                        remove_pr       = 1 - 0.976083968, 
                        lambda          = 0, 
                        res_death_type  = 1,
                        observe_type    = 2, #but only happens once a year
                        res_move_obs    = FALSE, 
                        max_ages        = 24, 
                        RESOURCE_ini    = 50, 
                        culling         = FALSE,
                        land_ownership  = TRUE,
                        age_repr        = 1,
                        land_dim_1      = DIM_1,
                        land_dim_2      = DIM_2,
                        manage_target   = 2500,
                        get_res         = 'Full');
old_obs  <- sim_old[["observation_array"]];

# BD: Add crop types on layer 1 of the landscape with the crop_type function
sim_old[["LAND"]] <- crop_type(land = sim_old[["LAND"]]);

# BD: Start out just recording 3 variables as output, but can increase
# BD: I've changed the notation to double-brackets, which is technically
#     a bit more secure. It works if you edit columns of the resource_array
#     directly.
sim_sum_1 <- matrix(data = NA, nrow = TMAX, ncol = 3);
for(time_step in 1:TMAX){
  sim_new                 <- gmse_apply(old_list = sim_old, get_res = 'Full');
  sim_sum_1[time_step, 1] <- time_step;
  sim_sum_1[time_step, 2] <- sim_new$basic_output$resource_results[1];
  sim_sum_1[time_step, 3] <- sim_new$basic_output$observation_results[1];
  next_time               <- time_step + 1; # What happens next time step?
  if(next_time %in% JAN){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.976083968; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% FEB){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.94824969; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% MAR){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.9023544; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% APR){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.811569975; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% MAY){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.805694059; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% JUN){
    sim_new[["resource_array"]][, 7]  <- 100;       # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.794513672; # Death probability
    sim_new[["resource_array"]][, 10] <- 1.18;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% JUL){
    sim_new[["resource_array"]][, 7]  <- 100;       # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.767700722; # Death probability
    sim_new[["resource_array"]][, 10] <- 1.18;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    temp_res                          <- sim_new[["resource_array"]];
    sim_new[["resource_array"]]       <- cell_K(res   = temp_res, DIM_1 = DIM_1, 
                                                DIM_2 = DIM_2);
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                      res  = temp_res);
  }
  if(next_time %in% AUG){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.842305069; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    old_obs                           <- sim_new[["observation_array"]];
    temp_res                          <- sim_new[["resource_array"]];
    sim_new[["resource_array"]]       <- cell_K(res   = temp_res, DIM_1 = DIM_1, 
                                                DIM_2 = DIM_2);
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_2(land = temp_land,  #function 2 because birth rate changes for type 3
                                                        res  = temp_res)
  }
  if(next_time %in% SEP){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.877005333; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    

    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  
  if(next_time %in% OCT){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.891729601; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
    
    
  }
  if(next_time %in% NOV){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.939895976; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  if(next_time %in% DEC){
    sim_new[["resource_array"]][, 7]  <- 0;         # Movement distance
    sim_new[["resource_array"]][, 9]  <- 1 - 0.956833861; # Death probability
    sim_new[["resource_array"]][, 10] <- 0;         # Birth probability
    sim_new[["observation_array"]]    <- old_obs;   # Use old observations
    
    temp_res                          <- sim_new[["resource_array"]];
    temp_land                         <- sim_new[["LAND"]];
    sim_new[["resource_array"]]       <- crop_hamster_1(land = temp_land, 
                                                        res  = temp_res);
  }
  sim_old <- sim_new; # BD: This should always go at the end
  
  print(sim_sum_1[time_step,]); # BD: Just to see the simulation progress
}

colnames(sim_sum_1) <- c("Time", "Pop_size", "Pop_est");
print(sim_sum_1) 




