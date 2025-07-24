
################ R_Script2 ###################

#Load libraries
library(here)
library(readxl)
library(httr)
library(xml2)
library(readr)
library(dplyr)
library(tidyverse)
library(purrr)
library(stringr)
library(fuzzyjoin)
library(tibble)

##Duplication Removal

#Two-step process was used for deduplication: String-matching algorithm using `synthesisr` package; Rayyan duplication removal tool (https://www.rayyan.ai/)

#String-matching algorithm using synthesisr
#load and process RIS file
ris_dat<- read_refs(here("search/Pest_ETH.ris"))
#1539 articles
ris_dat$title2 <- str_replace_all(ris_dat$title,"[:punct:]","") %>% 
  str_replace_all(.,"[ ]+", " ") %>% tolower()
ris_dat <- distinct(ris_dat, title2, .keep_all = TRUE)
#dim(ris_dat)

#remove partial matches in titles
duplicates_string <- find_duplicates(ris_dat$title2, method = "string_osa", to_lower = TRUE, rm_punctuation = TRUE, threshold = 7)
#manual_checks <- review_duplicates(ris_dat$title, duplicates_string)

#extract unique references and drop title2 & n_duplicates
ris_dat1 <- extract_unique_references(ris_dat, duplicates_string) %>%
  dplyr::select(-title2, -n_duplicates)
#1303 after duplicate removal!

#save as a .bib file, upload to Zotero, export as .ris file, import to Rayyan.
#write_refs(ris_dat1, format = "bib", file = "search/Pest_ETH_dedupl.bib")


## Supplementary Databases

#We need two ancillary/supplemental datasets in addition to our main dataset to enrich our analysis and answer our research questions (RQ), including pesticides' Maximum Residue Limit (MRL) and pesticides' Toxicity Relevance Value (TRV; chronic ADI/RfD).

#For MRL dataset, we used the publicly available EU MRL and USDA databases. For TRV dataset, we used EU pesticide properties (EFSA) and US EPA IRIS databases.  In both datasets, we preferred the lowest values whenever multiple values retrieved.

#Load the main dataset and select distinct pesticides
raw <- read_excel(here("Data/source_data.xlsx"), sheet = "pest_raw") %>%
  mutate(pest_met = str_to_lower(str_squish(pest_met)), 
         food = str_to_lower(str_squish(food)))

pest_mrl<- raw %>% distinct(pest_met, food) %>% #pesticide metabolites
  mutate(pest_met = str_to_lower(str_squish(pest_met)), 
         food = str_to_lower(str_squish(food))) %>%
  rename(pesticide = pest_met)

pest_trv<- raw %>% distinct(pest_met) %>% #pesticide metabolites
  mutate(pest_met = str_to_lower(str_squish(pest_met))) %>%
  rename(pesticide = pest_met)

## MRL dataset processing

#EU MRL database: We downloaded the recent complete EU MRL database in XML format from EU Pesticides Database (v3.3; https://ec.europa.eu/food/plant/pesticides/eu-pesticides-database/start/screen/mrls/download) and saved a folder "./pesticide_xmls" We then parsed each files using xml2 R package and extracted pesticide names, food product categories, and the corresponding MRL values.

#Convert xml files to csv
folder <- "Data/pesticide_xmls"
files  <- file.path(folder, paste0("Publication", 1:8, ".xml"))

#Function to parse one publication
parse_pub <- function(file) {
  doc  <- read_xml(file)
  subs <- xml_find_all(doc, "//Substances")
  # for each Substance node, extract metadata + products
  df_list <- lapply(subs, function(sub_node) {
    sub_name    <- xml_text(xml_find_first(sub_node, "./Name"))
    pest_resid  <- xml_text(xml_find_first(sub_node, "./Pest_res_id"))
    entry_force <- xml_text(xml_find_first(sub_node, "./entry_force"))
    directive   <- xml_text(xml_find_first(sub_node, "./directive"))
    
    products <- xml_find_all(sub_node, "./Product")
    tibble(Publication = basename(file),
           Substance = sub_name,
           Pest_res_id = pest_resid,
           Entry_force = entry_force,
           Directive = directive,
           Crop = xml_text(xml_find_all(products, "./Product_name")),
           Crop_code = xml_text(xml_find_all(products, "./Product_code")),
           MRL = xml_text(xml_find_all(products, "./MRL")),
           ApplicationDate = xml_text(xml_find_all(products, "./ApplicationDate")))})
  bind_rows(df_list)}

#Parse all and combine
mrl_eu <- bind_rows(lapply(files, parse_pub)) #~251990 obs.
#write.csv(mrl_eu, "Data/inputs/mrl_eu.csv", row.names = FALSE)

#Process csv file
mrl_eu <- mrl_eu %>%
  rename(pesticide = Substance, food = Crop, mrl = MRL) %>%
  mutate(pesticide = str_trim(pesticide),
         food = str_trim(food),
         pesticide = str_to_lower(str_squish(pesticide)),
         food = str_to_lower(str_squish(food)),
         mrl = as.numeric(str_remove(mrl, "\\*"))) %>%
  filter(!is.na(mrl))

#USDA MRL database:
mrl_usda <- read_excel(here("Data/inputs/mrl_usda.xlsx"), sheet = "Pesticide-Raw") %>%
  rename(pesticide = `Index Pesticide`, 
         food = `Index Commodity`, mrl = `MRL (ppm)`) %>%
  mutate(pesticide = str_to_lower(str_squish(pesticide)),
         food = str_to_lower(str_squish(food)),
         mrl = as.numeric(mrl))

## Match pesticide names
#Fuzzy Matching ~ EU MRL
mrl.pest.eu <- pest_mrl %>% distinct(pesticide) %>%
  stringdist_left_join(mrl_eu %>% distinct(pesticide),
                       by = "pesticide",
                       method = "jw", max_dist = 0.15, distance_col = "dist") %>%
  group_by(pesticide.x) %>%
  slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(pesticide.raw = pesticide.x, pesticide.mrl = pesticide.y)

#Assess Matching Results
mrl.pest.eu %>% mutate(matched = !is.na(pesticide.mrl)) %>%
  summarise(total_pests = n(), matched = sum(matched),
            pct_matched    = matched / total_pests * 100)
#Only 32/88 pesticides were matched.

#Fuzzy Matching ~ USDA MRL
mrl.pest.usda <- pest_mrl %>% distinct(pesticide) %>%
  stringdist_left_join(mrl_usda %>% distinct(pesticide),
                       by = c ("pesticide" = "pesticide"),
                       method = "jw", max_dist = 0.15, distance_col = "dist") %>%
  group_by(pesticide.x) %>% 
  slice_min(dist, n = 1, with_ties = FALSE) %>% ungroup() %>% 
  rename(pesticide.raw1 = pesticide.x, pesticide.mrl1 = pesticide.y)

#Assess Matching Results
mrl.pest.usda %>% mutate(matched = !is.na(pesticide.mrl1)) %>%
  summarise(total_pests = n(), matched = sum(matched),
            pct_matched    = matched / total_pests * 100)
#only 41/88 pesticides were matched.

#Note: The majority of the pesticides name across both dataset doesn't match. Thus, it better to proceed with manual matching.

#manually matched names and add to our main dataset
#pest_manual <- read.csv(here("Data/mrl_pest_name.csv"), na.strings = "") %>%
#  mutate(pesticide.mrl = str_to_lower(pesticide.mrl),
#         pesticide.mrl1 = str_to_lower(pesticide.mrl1))

#pest_mrl <- pest_mrl %>%
#  left_join(pest_manual, by = c("pesticide" = "pesticide.raw"))

#Fuzzy Match Food Subgroups ~ EU MRL
mrl.food.eu <- pest_mrl %>% distinct(food) %>%
  stringdist_left_join(mrl_eu %>% distinct(food),
                       by = "food",
                       method = "jw", max_dist = 0.15, distance_col = "dist") %>%
  group_by(food.x) %>%
  slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(food.raw = food.x, food.mrl = food.y)

#Assess Matching Results
mrl.food.eu %>% mutate(matched = !is.na(food.mrl)) %>%
  summarise(total_pests = n(), matched = sum(matched),
            pct_matched    = matched / total_pests * 100)
#only 8/18 foods were matched. 

#Fuzzy Match Food Subgroups ~ USDA MRL
mrl.food.usda <- pest_mrl %>% distinct(food) %>%
  stringdist_left_join(mrl_usda %>% distinct(food),
                       by = "food",
                       method = "jw", max_dist = 0.15, distance_col = "dist") %>%
  group_by(food.x) %>%
  slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(food.raw1 = food.x, food.mrl1 = food.y)

#Assess Matching Results
mrl.food.usda %>% mutate(matched = !is.na(food.mrl1)) %>%
  summarise(total_pests = n(), matched = sum(matched),
            pct_matched    = matched / total_pests * 100)
#only 8/18 foods were matched. 

#Let's do the matching for all the available foods manually.
food_manual <- tibble::tibble(
  food.raw = c("cabbage","corn","drinking water","fish","honey","khat","meat", "milk","millet","onion","pineapple","potato","rice","sorghum", "swiss chard","tea","tomato","wheat"),
  
  food.mrl = c("head cabbages","maize/corn", NA,NA, "honey and other apiculture products (7)","(b) leaves and herbs","muscle","milk", "common millet/proso millet","onions","pineapples","(a) potatoes","rice","sorghum", "chards/beet leaves","teas","tomatoes","wheat"),
  
  food.mrl1 = c("cabbage","corn, grain", NA, "catfish, freshwater", "honey", "swiss chard","cattle, meat","milk", "millet, pearl, grain", "onion, bulb","pineapple","potato","rice","sorghum, grain", "swiss chard","tea, leaves","tomato","wheat, grain"))

#manually matched names and add to our main dataset
pest_mrl <- pest_mrl %>%
  left_join(food_manual, by = c("food" = "food.raw"))

## Final pesticide-food-MRL dataset

#Keep the row with the minimum MRL value for each pesticide-food group
mrl_eu.low <- mrl_eu %>%
  group_by(pesticide, food) %>%
  slice(which.min(as.numeric(mrl))) %>%
  ungroup()

mrl_usda.low <- mrl_usda %>%
  group_by(pesticide, food) %>%
  slice(which.min(as.numeric(mrl))) %>%
  ungroup()

pest_mrl1<- pest_mrl %>%
  left_join(mrl_eu.low %>% select(pesticide, food, mrl), 
            by = c("pesticide.mrl" = "pesticide", "food.mrl" = "food")) %>% 
  left_join(mrl_usda.low %>% select(pesticide, food, mrl), 
            by = c("pesticide.mrl1" = "pesticide", "food.mrl1" = "food")) %>%
  rename(mrl_eu = mrl.x, mrl_usda = mrl.y)

#For cases where EU and USDA MRL present, select the lowest MRL
pest_mrl1 <- pest_mrl1 %>%
  mutate(MRL = case_when(
    !is.na(mrl_eu) & !is.na(mrl_usda) ~ pmin(as.numeric(mrl_eu), as.numeric(mrl_usda), na.rm = TRUE),
      is.na(mrl_eu) & !is.na(mrl_usda) ~ as.numeric(mrl_usda),
      !is.na(mrl_eu) & is.na(mrl_usda) ~ as.numeric(mrl_eu),
      TRUE ~ NA_real_))

#Apply default MRL value (0.1 µg/L) for drinking water
pest_mrl1 <- pest_mrl1 %>%
  mutate(MRL = if_else(is.na(MRL) & str_detect(food, "drinking water"), 
                       0.1, MRL))

#Re-check final MRL matching
pest_mrl1 %>% mutate(matched = !is.na(MRL)) %>%
  summarise(total = n(), matched = sum(matched),
            pct_matched    = matched / total * 100)
#277/359 mrls were finally retrieved

#check manually for missing mrls
pest_mrl_check <- pest_mrl1 %>% distinct(pesticide, food, MRL) %>%
  filter(is.na(MRL))

#Re-check if available in Codex Pesticides Residues in Food Online Database

#save the dataset
pest_mrl_df <- pest_mrl1 %>% 
  select(pesticide, food, MRL)

write.csv(pest_mrl_df, "Data/inputs/pest_mrl_df.csv", row.names = FALSE)


## TRV dataset processing

#EU active pesticides database: We downloaded the recent complete EU MRL database in XML format from EU Pesticides Database (v3.3; https://ec.europa.eu/food/plant/pesticides/eu-pesticides-database/start/screen/mrls/download) and saved a folder "./pesticide_xmls" We then parsed each files using xml2 R package and extracted pesticide names, food product categories, and the corresponding MRL values.

#select the chronic TRV ~ ADI (mg/kg bw/d)
trv_eu <- read_excel(here("Data/inputs/trv_efsa.xlsx")) %>%
  rename(pesticide = Substance, ADI = `ADI (mg/kg bw/d)`) %>%
  mutate(pesticide = str_to_lower(str_squish(pesticide)),
         ADI = as.numeric(str_remove(ADI, "\\mg/kg bw/day"))) %>%
  dplyr::select(pesticide, ADI) 

#US EPA IRIS database:
#We download the complete IRIS database (last update:15/04/2025; https://www.epa.gov/system/files/documents/2025-04/iris_downloads_database_export_april2025.xlsx) and processed it as follow.

#Load and clean the data: TRV ~ RfD (mg/kg-day)
trv_iris_rfd <- read_excel(here("Data/inputs/trv_iris.xlsx"), sheet = "RfD Toxicity Values") %>%
  rename(pesticide = `CHEMICAL NAME`,
         RfD = `RFD VALUE`,
         quality = `OVERALL CONFIDENCE`) %>%
  dplyr::select(pesticide, RfD, quality) %>%
  mutate(pesticide = str_to_lower(str_squish(pesticide)),
         RfD = as.numeric(RfD),
         quality = case_when(quality == "Low/Medium" ~ "Medium-Low",
                             quality == "Medium/High" ~ "Medium-High",
                             TRUE ~ quality))

#If multiple RfD available for a single pesticide, select the one with better quality and if quality is the same, select the lowest value "most protective" approach ~ all irrespective of effect type!!

#order quality from worst to best quality
quality_order <- c("Low", "Medium-Low", "Medium", "Medium-High", "High")

#select based on quality and then lowest RfD
trv_iris_rfd <- trv_iris_rfd %>%
  arrange(pesticide,
          factor(quality, levels = quality_order, ordered = TRUE) %>% desc(),
          RfD) %>%
  distinct(pesticide, .keep_all = TRUE) %>%
  select(-quality)

#Load and clean the data: TRV ~ OSF (per mg/kg-day)
trv_iris_osf <- read_excel(here("inputs/iris_downloads_database_export_april2025.xlsx"), sheet = "WOE Toxicity Values") %>%
  rename(pesticide = `CHEMICAL NAME`,
         route = `CA ROUTE`,
         OSF = `SLOPEFACTOR_VALUE`) %>%
  select(pesticide, route, OSF) %>%
  mutate(pesticide = str_to_lower(str_squish(pesticide)),
         OSF = as.numeric(OSF)) %>%
  filter(route == "Oral") %>% select(-route)

## Final Matching ~ add TRV (ADI, RfD & OSF) to our dataset

#match pesticides name manually
pest_trv_man <- read_excel(here("inputs/pest_name_trv.xlsx")) %>% 
  mutate(pesticide.trv = str_to_lower(str_squish(pesticide.trv)), 
         pesticide.trv1 = str_to_lower(str_squish(pesticide.trv1)))

#pest_trv<- pest_trv %>%
#  left_join(pest_trv_man, by = c("pesticide" = "pesticide.raw"))

#add EFSA and IRIS TRV
pest_trv <- pest_trv_man %>%
  left_join(trv_eu, by = c("pesticide.trv" = "pesticide")) %>% 
  left_join(trv_iris_rfd, by = c("pesticide.trv1" = "pesticide")) %>% 
  left_join(trv_iris_osf, by = c("pesticide.trv1" = "pesticide"))

#assess matching results
sum(!is.na(pest_trv$ADI)) #retrieved for 39 pests
sum(!is.na(pest_trv$RfD)) #retrieved for 47 pests
sum(!is.na(pest_trv$OSF)) #retrieved for 23 pests

#merge ADI and RfD as chronic TRV by selecting the lowest value
pest_trv <- pest_trv %>%
  mutate(TRV = case_when(
    !is.na(RfD) & (RfD <= ADI | is.na(ADI)) ~ RfD, 
    !is.na(ADI) ~ ADI, TRUE ~ NA_real_),
    source1 = case_when(!is.na(RfD) & (RfD <= ADI | is.na(ADI)) ~ "EPA", 
                       !is.na(ADI) ~ "EFSA", TRUE ~ NA_character_),
    source2 = case_when(!is.na(OSF) ~ "EPA", TRUE ~ NA_character_)) %>%
  rename(pesticide = pesticide.raw) %>% 
  select(pesticide, TRV, source1, OSF, source2)

#find missing TRVs from other sources
pest_trv_m <- pest_trv %>% filter(is.na(OSF)) %>%
  select(pesticide)

#retrieved from PPDB: Pesticide Properties DataBase (https://sitem.herts.ac.uk/aeru/ppdb/en/index.htm)

pest_trv1 <- pest_trv %>%
  mutate(
    TRV = case_when(
      is.na(TRV) & pesticide == "chlorpyrifos" ~ 0.001,
      is.na(TRV) & pesticide == "dde" ~ 5.0e-04,
      is.na(TRV) & pesticide == "difenocozole" ~ 0.01,
      is.na(TRV) & pesticide == "fenthion" ~ 0.01, #ARfD
      is.na(TRV) & pesticide == "hexaconazole" ~ 0.005,
      is.na(TRV) & pesticide == "o,p'-ddd" ~ 5.0e-04,
      is.na(TRV) & pesticide == "o,p'-dde" ~ 5.0e-04,
      is.na(TRV) & pesticide == "p,p'-ddd" ~ 5.0e-04,
      is.na(TRV) & pesticide == "p,p'-dde" ~ 5.0e-04,
      is.na(TRV) & pesticide == "piperonyl-butoxide" ~ 0.2,
      is.na(TRV) & pesticide == "tebuconazole" ~ 0.03,
      TRUE ~ TRV),
    source1 = case_when(
      pesticide == "chlorpyrifos" & is.na(source1) ~ "PPDB",
      pesticide == "dde" & is.na(source1) ~ "EPA",
      pesticide == "difenocozole" & is.na(source1) ~ "PPDB",
      pesticide == "fenthion" & is.na(source1) ~ "PPDB",
      pesticide == "hexaconazole" & is.na(source1) ~ "PPDB",
      pesticide == "o,p'-ddd" & is.na(source1) ~ "EPA",
      pesticide == "o,p'-dde" & is.na(source1) ~ "EPA",
      pesticide == "p,p'-ddd" & is.na(source1) ~ "EPA",
      pesticide == "p,p'-dde" & is.na(source1) ~ "EPA",
      pesticide == "piperonyl-butoxide" & is.na(source1) ~ "PPDB",
      pesticide == "tebuconazole" & is.na(source1) ~ "PPDB",
      TRUE ~ source1)
  )

#write.csv(pest_trv1, "Data/inputs/pest_trv_df.csv", row.names = FALSE)

##NOTE: In cases where multiple TRV found for the same pesticide/metabolite, we selected the lowest TRV ("most protective" approach). Similarly, if single TRV is found for only parent pesticides/metabolites, it was used for all metabolites related to that pesticide. If multiple TRV is found for different metabolites and missing for related other metabolite, TRV of more related metabolite was used (vv. parent pesticides).

## Food consumption

#Ethiopian food consumption data was retrieved from the publicly available Ethiopia Socioeconomic Survey fourth round panel dataset (ESS4; 2018/2019) [https://doi.org/10.48529/k739-c548]. In comparison to previous round surveys, ESS4 has broader coverage and considered more representative at national level. It covered all regional states and two administrative cities (Addis Ababa and Dire Dawa) with a total of 565 enumeration areas, of which 316 are rural and 219 are urban. ESS4 included 7,527 households selected through multi-stage sampling procedures. After downloading the dataset, specific Food consumption data file (sect6a_hh_w4.dta) was selected, imported, cleaned and processed as below. The data file included household food consumption (quantity and value) in the last 7 days and source of foods consumed by the household from a subset list of food items. As per the eligibility criteria applied to EES4, we used a total of n=20932 unique survey results exactly matching foods in our dataset (n=17, except for drinking water) to calculate the overall national food consumption rate (mean, min, max, 99th percentile).

#citation:
#Central Statistical Agency of Ethiopia. Ethiopia Socioeconomic Survey (ESS4) 2018-2019. Public Use Dataset. Ref: ETH_2018_ESS_v03. Downloaded from[https://microdata.worldbank.org/index.php/catalog/3823/get-microdata] on [03/06/2025]. DOI [https://doi.org/10.48529/k739-c548]

#dataset cleaning and processing
#Our interest is the amount of food consumed daily (kg/day). Select food group, amount & unit; where foods should match our main dataset, consumption rates are presented and standard measurement units are employed.

fc_raw <- read.csv(here("Data/inputs/fc_raw.csv"),na.strings = "") %>%
  rename(food = item_cd, area = saq14, region = saq01, 
         amount = s6aq02a, unit = s6aq02b) %>%
  dplyr::select(food, area, region, amount, unit) %>%
  mutate(food = str_replace(food, "^\\d+\\.\\s*", ""),
         area = str_replace(area, "^\\d+\\.\\s*", ""),
         region = str_replace(region, "^\\d+\\.\\s*", ""),
         unit = str_replace(unit, "^\\d+\\.\\s*", "")) %>%
  filter(!is.na(amount)) #drop missing data

#fc_raw_type <- fc_raw %>% distinct(food, food2)
#fc_raw_unit <- fc_raw %>% distinct(unit1, unit2) #use kilogram/gram

#Match with and select foods that are only presented in our dataset. If further food subgroups were presented (disaggregated food groups), aggregate to match our dataset, e.g., beef, goat & mutton meat can be aggregated as meat.

fc_raw <- fc_raw %>%
  mutate(food = case_when(
      food == "Goat & mutton meat" ~ "Meat",
      food == "Beef" ~ "Meat",
      food == "Maize" ~ "Corn",
      food == "Wheat (Incl. Flour factory product)" ~ "Wheat",
      food == "Honey, natural" ~ "Honey",
      food == "Chat / Kat" ~ "Khat",
      food == "Mango" ~ "Pineapple",  #assumed based on local knowledge
      food == "kale, cabbage, Pumpikn Leaf, Lettuce, spinach" ~ "Cabbage",
      food == "kale, cabbage, Pumpikn Leaf, Lettuce, spinach" ~ "swiss chard",
      TRUE ~ food),
      food = str_to_lower(str_squish(food)))

fc_manual <- raw %>% distinct(food) %>%
  stringdist_left_join(fc_raw %>% distinct(food), by = "food",
            method = "jw", max_dist = 0.1, distance_col = "dist") %>%
  group_by(food.x) %>% slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>% rename(food.raw = food.x, food.nw = food.y)

#select food groups in our data
fc_raw1 <- fc_raw %>%
  filter(food %in% fc_manual$food.nw) %>%
  filter(unit == "Gram" | unit == "Kilogram")

# Convert amount in Gram to Kilogram
fc_raw1 <- fc_raw1 %>%
  mutate(amount = as.numeric(amount),
         amount = case_when(unit == "Gram" ~ amount / 1000,
                            unit == "Kilogram" ~ amount, 
                            TRUE ~ NA_real_))

#Final food consumption rate (CRi) calculations
fc_raw_f <- fc_raw1 %>% group_by(food) %>%
  summarise(CRi = round (mean(amount),2),  #average
            CRi_sd = round (sd(amount),2),
            CRi_min = round(min(amount),2),
            CRi_max = round(max(amount),2),
            CRi_P99 = round(quantile(amount, 0.99),2))

# Add missing groups (drinking water and swiss chard)
fc_raw_f <- bind_rows(
  fc_raw_f,
  data.frame(food = "drinking water", CRi = NA_real_),
  data.frame(food = "swiss chard", CRi = NA_real_))

fc_raw_f <- fc_raw_f %>%
  mutate(
    CRi = case_when(
      food == "drinking water" ~ 1.95,  #standard water consumption
      food == "swiss chard" ~ 1.79,     #~kale, cabbage, Lettuce, spinach etc.
      TRUE ~ CRi),
    CRi_sd = case_when(
      food == "drinking water" ~ 0.64,  
      food == "swiss chard" ~ 1.94,     
      TRUE ~ CRi))

#write.csv(fc_raw_f, "Data/inputs/fc_processed.csv", row.names = FALSE)


#Identify best fit for Ci (n=2271) & CRi (n=20932) distribution

res_dat1 <- pest_df_comp %>%
  dplyr::select(pest, food, log_Mean, log_SD) %>%
  mutate(pesticide = pest, Mean = exp(log_Mean), SD = exp(log_SD)) %>%
  dplyr::select(-log_Mean, -log_SD)

fc_dat1<- fc_raw1 %>% dplyr::select(food, amount)

#Select fit variable
residue_values <- res_dat1$Mean
consumption_values <- fc_dat1$amount

#Fit distributions
fit_ln <- fitdist(residue_values, "lnorm")
fit_gamma <- fitdist(residue_values, "gamma")
fit_weibull <- fitdist(residue_values, "weibull")
fit_exp <- fitdist(residue_values, "exp")

fit_ln2 <- fitdist(consumption_values, "lnorm")
fit_gamma2 <- fitdist(consumption_values, "gamma")
fit_weibull2 <- fitdist(consumption_values, "weibull")
fit_exp2 <- fitdist(consumption_values, "exp")

#Compare fits
gof_residue <- gofstat(
  list(fit_ln, fit_gamma, fit_weibull, fit_exp),
  fitnames = c("Lognormal", "Gamma", "Weibull", "Exponential"))

gof_consumption <- gofstat(
  list(fit_ln2, fit_gamma2, fit_weibull2, fit_exp2),
  fitnames = c("Lognormal", "Gamma", "Weibull", "Exponential"))

print(gof_residue)
print(gof_consumption)
#NOTE: Lognormal distribution was found as best best-fit for both variables (~Best AIC/BIC and GOF overall)! Now, we use lognormal distribution for simulations of point estimates of pooled pesticide residue concentration (Mean [95% CI]) and Consumption Rates (mean, min, max/99 percentile) to calculate risks according for each pesticide-food combination.

##____________________//END//________________________##
