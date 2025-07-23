# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
#------------------------------#
# STEP 1: Load and merge WGS batch and mtDNA-CN data
#------------------------------#
# mtDNA-CN summary table across all sequenced UKB participants
ukbdata <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/mitochondria/20240423_mtDNA_ukbb_batch1_batch2_summary.csv") # mtDNA-CN (n = 485,451)
# Sequencing batch assignment
batch <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/mitochondria/20240423_UKBB_eid_sequencing_batch.tsv") # Batch info: Batch 1 (2021 release) = 200,017, Batch 2 (2023 release) = 285,434
# Merge batch info
ukbdata <- ukbdata %>% 
  full_join(batch,by=c("sample"="eid")) # mtDNA-CN (n = 485,451)
# Exclude withdrawn participants
withdrawn <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/withdrawal_notifications/w36610_2023-04-25.csv") # withdrawn participants (n = 162)
ukbdata <- ukbdata %>% 
  filter(!(sample %in% withdrawn$V1)) # mtDNA-CN (n = 485,443)
#------------------------------#
# mtDNA-CN Raw data
mt_raw <- ukbdata %>% dplyr::select(sample,mt) #mtDNA-CN Raw (n = 485,443)
#------------------------------#
# Select relevant fields & define PD status
ukbdata <- ukbdata %>% 
  dplyr::select(sample,mt,genetic_sex,Blood_Collection_Age,logmt,batch,PC1,PC2,PC3,PC4,PC5,Parkinsons) %>% 
  mutate(PD_diagnosis = ifelse(Parkinsons=="Y","Case","Control")) %>% 
  dplyr::select(-Parkinsons) %>% 
  na.omit()
# Rename for clarity
colnames(ukbdata) <- c("sample", "mt", "sex", "age","logmt", "Batch", "PC1","PC2", "PC3","PC4", "PC5", "diagnosis") # mtDNA-CN (n = 432,901)
#------------------------------#
# STEP 2: Load and clean blood count data
#------------------------------#
# Load blood measurements (baseline measurement instance 0 only)
blood_raw <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/extracts/20240118_673741_ukb_blood.csv.gz") # blood measurements (n = 502,366)
blood <- blood_raw %>%
  dplyr::select(eid, matches("0-0.0")) %>%
  dplyr::select(eid, where(is.numeric))
# Pivot long
blood_long <- blood %>%
  pivot_longer(!eid, names_to = "FieldIDNew", values_to = "Value")

# Map UKB field names
dictionary <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/extracts/20240219_Data_Dictionary_Showcase.csv") %>%
  mutate(FieldIDNew = paste0(FieldID, "-0.0")) %>%
  dplyr::select(FieldIDNew, Field)
# Join and pivot wide to get blood measurement names
blood_named <- left_join(blood_long, dictionary, by = "FieldIDNew") %>%
  select(eid, Field, Value) %>%
  pivot_wider(names_from = Field, values_from = Value)

# Merge blood measurement to mtDNA data
ukbdata <- left_join(ukbdata, blood_named, by = c("sample" = "eid")) # n = 432,901

#------------------------------#
# STEP 3: Remove unwanted blood features
#------------------------------#
# Remove blood traits with >1% missing
miss_rate <- colSums(is.na(ukbdata)) / nrow(ukbdata) * 100
ukbdata <- ukbdata[, names(miss_rate[miss_rate <= 1]), with = FALSE]

# Exclude low-abundance cell types and RBC traits (no mtDNA)
exclude_vars <- c(
  "Basophill percentage", "Nucleated red blood cell percentage",
  "Basophill count", "Nucleated red blood cell count",
  "Red blood cell (erythrocyte) count", "Haemoglobin concentration", "Haematocrit percentage",
  "Mean corpuscular volume", "Mean corpuscular haemoglobin", "Mean corpuscular haemoglobin concentration",
  "Red blood cell (erythrocyte) distribution width"
)

ukbdata <- ukbdata %>% select(-any_of(exclude_vars)) # n = 432,901
#------------------------------#
# STEP 4: Select blood traits for model adjustment
#------------------------------#
# Final selected blood cell traits
final_traits <- c(
  "White blood cell (leukocyte) count", "Platelet count", "Platelet crit",
  "Mean platelet (thrombocyte) volume", "Platelet distribution width",
  "Lymphocyte percentage", "Monocyte percentage",
  "Neutrophill percentage", "Eosinophill percentage"
) # Use blood cell percentages when available, as they better fit the model. For traits without percentage values, counts are used instead.

# Subset data for model
blood_model_data <- ukbdata %>%
  dplyr::select(sample, mt, sex, age, logmt, Batch, PC1, PC2, PC3, PC4, PC5, diagnosis, all_of(final_traits)) # n = 432,901

#------------------------------#
# STEP 5: Exclude related individuals & missing
#------------------------------#
# Remove 3rd-degree relatives (kinship > 0.0625)
related <- fread("/vast/projects/bahlo_ukbiobank/app36610_HDA/mitochondria/20240426_mwvc_pd_todrop.txt_MWVC.txt") # n = 63,771
blood_model_data <- blood_model_data %>% filter(!(sample %in% related$V1)) %>% na.omit() #Final sample size n = 367,322
#------------------------------#
# STEP 6: Stepwise model selection
#------------------------------#
# Build linear model predicting logmt ~ blood composition
model_input <- blood_model_data %>%
  select(logmt, all_of(final_traits))
# Fit initial linear model
m0 <- lm(logmt ~ ., data = model_input)
# Perform stepwise selection to identify best-fitting model
m1 <- step(m0, direction = "both")  # Stepwise AIC model selection
summary(m1)  # Adjusted R² ~19%
#------------------------------#
# STEP 7: Generate adjusted mtDNA-CN (residuals)
#------------------------------#
# Residuals = adjusted mtDNA-CN (removes blood composition effects)
blood_model_data <- blood_model_data %>%
  mutate(mt.adj = round(residuals(m1), 4))
#------------------------------#
# mtDNA-CN adjusted data
mt_adjust <- blood_model_data %>% dplyr::select(sample,mt.adj) #mtDNA-CN adjusted data (n = 367,322)
#------------------------------#
# Save final data
mtDNA_CN <- mt_raw %>% 
  full_join(mt_adjust) # mtDNA-CN Raw (n = 485,443); mtDNA-CN adjusted data (n = 367,322)
colnames(mtDNA_CN) <- c("eid","mtDNA_CN_raw","mtDNA_CN_adjusted")
# After the paper was published, a new withdrawal list was released (Dec 2024)
# Update the dataset by excluding these newly withdrawn participants
new_withdraw <- read.csv("/vast/projects/bahlo_ukbiobank/app36610_HDA/withdrawal_notifications/w36610_20241217.csv", header = FALSE) # n = 398
mtDNA_CN <- mtDNA_CN %>%
  filter(!(eid %in% new_withdraw$V1)) # mtDNA-CN Raw (n = 485,270); mtDNA-CN adjusted data (n = 367,206)
#write.csv(mtDNA_CN, "/vast/projects/bahlo_ukbiobank/app36610_HDA/mitochondria/UKBB_mtDNA-CN_Return/mtDNA_CN_UKBB.csv", row.names = FALSE)

