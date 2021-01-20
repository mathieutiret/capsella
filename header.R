## loading (libraries, helper functions). -----------------
source( "~/work/general/src/R/include.R" )

## global variables. --------------------------------------
.dataset = c( "uppsala", "guangzhou", "toronto" )
.cluster = c( "ASI", "EUR", "ME" )
.ccl     = c( "#3dcd64", "red", "#7C40F5" ) %>% name( .cluster )
.nblock  = 4
.pheno   = c( "GP", "BP", "FP", "SP" )

## global dataset. ----------------------------------------
# phenotype data. -------------------------------
.data = 
  # loading the original dataset. ---------------
  read_tibble( here( "data", "common_garden.csv" ) ) %>%
  
  # correcting encoding errors. -----------------
  # standardising sample IDs. 
  mutate( sample_name = gsub( "_", "", sample_name ) ) %>% 

  # recoding. -----------------------------------
  # missing values.
  # @dev : automatic coercion with %as.numeric. 
  mutate_at(
    vars( matches("Mean_*"), height_FT, transfer_time ),
    ~ suppressWarnings( num(.) )
  ) %>% 
  
  # genetic clusters.
  mutate( 
    cluster_adm = 
      recode( cluster_admixture, 
              blue = "ME", green = "ASI", red = "EUR" ),
    cluster_mol = 
      recode( cluster_admixture_MOLECOL,
              blue = "ME", green = "ASI", red = "EUR" )
  ) %>% 
  
  mutate_at( vars(AccessionID), chr ) %>% 
  left_join(
    read_tibble( here( "data", "safe_cluster.csv" ) ) %>% 
      select( AccessionID, cluster_new )              %>% 
      mutate_at( vars(AccessionID), chr ),
    by = "AccessionID"
  ) %>%
  
  # countries & regions. 
  mutate(
    country = recode(
      country,
      Russian            = "Russia",
      Czech_Republic     = "Czech" ,
      Bosnia_Hercegovina = "Bosnia"
    ),
    region =
      case_when(
        country == "United_States" ~ "usa",
        country %in% c( "China", "Taiwan" ) ~ "asia",
        country %in% 
          c( "Algeria", "Jordan", "Israel", "Syria", "Turkey" ) ~ 
          "middle east",
        T ~ "europe"
      )
  ) %>%
  
  # block positions.
  mutate(
    row = str_extract( case_ID, "[A-Z]" ),
    col = str_extract( case_ID, "[0-9]+" )
  ) %>% 
  
  # adding variables. ---------------------------
  mutate(
    fitness = 
      ifelse(
        dataset == "toronto" | height_main_inflorescence < 100,
        number_fruits_branches_10cm,
        number_fruits_branches_10cm / 100 * 
          nb_primary_branchs * height_main_inflorescence  
      ),
    
    GP = germination_time,
    BP = bolting_time - germination_time,
    FP = flowering_day - bolting_time,
    SP = senescence_date - flowering_day
  ) %>% 
  
  # selecting & renaming variables. -------------
  select(
    # IDs.
    accessionID = AccessionID, sampleID = sample_name, 
    dataset, country, region, site = pop, latitude, longitude,
    replicateID, blockID = box_ID, 
    row, col,
    fit_recorderID = recoderID_fitness, 
    TT_recorderID  = transferID, 
    GT_recorderID  = recorderID_germination, 
    BT_recorderID  = recorderID_bolding, 
    FT_recorderID  = recorderID_floweringT, 
    cluster_adm, cluster_mol, cluster_new, 
    
    # matches("bio*"),
    
    # phenology.
    # @dev : not loading `FT` and `BT`, they are problematic.
    TT = transfer_time, broken_primary_branch,
    GT = germination_time, GT_temp = Mean_temp_GT, GT_hum = Mean_hum_GT,
    BT = bolting_time, BT_temp = Mean_temp_BT, BT_hum = Mean_hum_BT,
    FT = flowering_day, FT_temp = Mean_temp_FT, FT_hum = Mean_hum_FT,
    nb_flowers_FT, height_FT, survival_to_flowering,
    ST = senescence_date, ST_temp = Mean_temp_SENESCENCE, 
    ST_hum = Mean_hum_SENESCENCE,
    nb_l_ros = nb_leaves_rosette, diam_ros = max_rosette_diameter, 
    GP, BP, FP, SP,
    
    # phenotypes.
    nb_pbr = nb_primary_branchs, nb_sinf = nb_secondary_inflo, 
    height = height_main_inflorescence,
    nb_f_10cm = number_fruits_branches_10cm, 
    f_width = fruit_width_mean_mm, f_height = fruit_height_mean_mm, 
    nb_seed = mean_number_seed, 
    env_pls = env_pls_upp_guang_tor,
    fitness
  ) 


# environmental data. ---------------------------
env.data = 
  # loading the original dataset. ---------------
  read_tibble( here( "data", "environment.csv" ) ) %>% 
  
  # correcting variables. -----------------------
  # day length.
  mutate(
    daylength_pptions_on24hours =
      case_when(
        site == "uppsala" ~ 
          daylength(59.85882, 119 + day_expe),
        site == "guangzhou" ~ 
          daylength(23.12911, (334 + day_expe) %% 365 + 1),
        site == "toronto" ~ 
          daylength(43.653226, 151 + day_expe )
      )
  ) %>% 

  # humidity in Toronto. 
  mutate(
    Mean_hum = 
    ifelse(
      site == "toronto",
      read_tibble( here("data", "environment_toronto.csv") ) %>% 
        pull(HumidityAvg),
      Mean_hum
    )
  ) %>% 
  
  # renaming variables. -------------------------
  rename(
    temperature = Mean_temp, 
    humidity    = Mean_hum,
    daylength   = daylength_pptions_on24hours
  ) %>% 
  
  # adding months. ------------------------------
  mutate(
    month = case_when(
      site == "uppsala" ~ 
        as.Date( day_expe, origin = "2014/05/01" ) %>% 
        str_extract_all( "\\d+", simplify = T ) %>% .[,2] %>% num,
      site == "guangzhou" ~ 
        as.Date( day_expe, origin = "2014/11/01" ) %>% 
        str_extract_all( "\\d+", simplify = T ) %>% .[,2] %>% num,
      site == "toronto" ~ 
        as.Date( day_expe, origin = "2014/06/02" ) %>% 
        str_extract_all( "\\d+", simplify = T ) %>% .[,2] %>% num
    )
  )


# genetic data. ---------------------------------
gen.data = 
  # loading the original data. ------------------
  read_tibble( here( "data", "gbs.csv" ) ) %>% 
  
  # selecting columns. --------------------------
  select( accessionID = AccessionID, matches("^loc") ) %>% 
  
  # tidying. ------------------------------------
  # @dev : faster than automatic coercion. 
  pivot_longer(
    cols = !accessionID,
    names_to = "locus", values_to = "genotype",
    values_transform = list( genotype = chr )
  ) %>% 

  # recoding missing values. --------------------
  mutate(
    accessionID = chr(accessionID),
    genotype = replace( genotype, genotype == "-", NA )
  ) %>% 

  # joining genetic clusters. -------------------
  left_join(
    .data %>% distinct( accessionID, cluster_new ),
    by = "accessionID"
  )


## global functions. --------------------------------------
balance_block = function( data, n = .nblock ) 
{
  # computing permutations. ---------------------
  block.perm = combn( data %>% distinct(blockID) %>% pull, n )
  
  # computing the largest permutation. ----------
  index = 
    block.perm %>% 
    apply( 2, 
           function(x)
             data                          %>% 
             group_by( accessionID )       %>% 
             filter( all(x %in% blockID) ) %>% 
             ungroup()                     %>% 
             distinct( accessionID )       %>% 
             count()                       %>% 
             pull()
           ) %>% 
    which.max()
  
  # filtering. ----------------------------------
  data %>% 
    # filtering for duplicates.
    # @dev : keeping only the first one. 
    group_by( accessionID, blockID ) %>% 
    slice(1)                         %>%
    ungroup()                        %>% 
    
    # filtering for the best combination.
    group_by( accessionID )          %>% 
    filter( 
      all( block.perm[,index] %in% blockID ),
      blockID %in% block.perm[,index] 
      )                              %>% 
    ungroup()
}

