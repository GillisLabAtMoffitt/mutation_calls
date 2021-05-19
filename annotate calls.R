# import library
library(tidyverse)


#######################################################################################  I  ### Load data----
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CIBMTR")

non_annotated_calls <-
  read_csv(paste0(path, "/05172021/CIBMTR_No Annotation Filtering_cleaned vcf.csv"))
extra_annotations <- 
  read.delim(paste0(path, "/05172021/MoffitCalls.annotated.txt"))

# Clean extra_annotations to create var for join
extra_annotations <- extra_annotations %>% 
  mutate(CHROM = str_match(Variant, "chr(.*?):")[,2]) %>% 
  mutate(POS = as.numeric(str_match(Variant, ":([[:digit:]]*)")[,2])) %>% 
  mutate(REF = str_match(Variant, ":[[:digit:]]*(.*)>")[,2]) %>% 
  mutate(ALT = str_match(Variant, ">(.*)")[,2])

# write_csv(extra_annotations, paste0(path, "/05172021/MoffittCalls.annotated separated " , format(Sys.time(),"%m%d%Y"), ".csv"))
# Submit this file to Franklin

# Prep rescued mutations from MoffitCalls.annotated analyzed by Franklin
rescue_mutations <-
  read.delim("/Users/colinccm/Documents/GitHub/Gillis/no vpn data/CIBMTR/05172021/export_1621355280244.tsv") %>% 
  filter(str_detect(Submission.Date, "2021-05-18")) %>% 
  # select("Chrom", "Position", "Ref", "Alt", "Gene", "HgvsC", "HgvsP", "Effect", "Region", "Genoox.Classification") %>% 
  mutate(HgvsP1 = str_replace_all(HgvsP, c(
    "Ala" = "A" ,
    "Arg" = "R" ,
    "Asn" = "N" ,
    "Asp" = "D" ,
    "Cys" = "C" ,
    "Glu" = "E" ,
    "Gln" = "Q" ,
    "Gly" = "G" ,
    "His" = "H" ,
    "Ile" = "I" ,
    "Leu" = "L" ,
    "Lys" = "K" ,
    "Met" = "M" ,
    "Phe" = "F" ,
    "Pro" = "P" ,
    "Ser" = "S" ,
    "Thr" = "T" ,
    "Trp" = "W" ,
    "Tyr" = "Y" ,
    "Val" = "V" 
  ))) %>% 
  mutate(HgvsP1 = str_remove(HgvsP1, "p.")) %>% 
  # select("Chrom", "Position", "Ref", "Alt", "Gene", "HgvsC", "HgvsP", "Effect", "Region", "Genoox.Classification", "HgvsP1") %>% 
  filter(
    (str_detect(Gene, "DNMT3A") & str_detect(HgvsP1, "R882C")) | 
      (str_detect(Gene, "DNMT3A") & str_detect(HgvsP1, "R882H")) | 
      (str_detect(Gene, "DNMT3A") & str_detect(HgvsP1, "P904L")) | 
      (str_detect(Gene, "GNB1") & str_detect(HgvsP1, "K57E")) | 
      (str_detect(Gene, "JAK2") & str_detect(HgvsP1, "V617F")) |
      (str_detect(Gene, "SF3B1") & str_detect(HgvsP1, "K666N")) | 
      (str_detect(Gene, "SF3B1") & str_detect(HgvsP1, "K700E")) | 
      (str_detect(Gene, "SFRS2") & str_detect(HgvsP1, "P95L")) | 
      (str_detect(Gene, "U2AF1") & str_detect(HgvsP1, "S34F")) | 
      (str_detect(Gene, "U2AF1") & str_detect(HgvsP1, "Q157P")) | # U2AF1L5, U2AF1???
      (str_detect(Gene, "U2AF1") & str_detect(HgvsP1, "Q157R"))
  ) %>% 
  mutate(Chrom = str_remove(Chrom, "chr")) #%>% 
  # mutate(merging_rescue = paste0(Chrom, Ref, Position, Alt))


# Join with non_annotated_calls with extra_annotations
non_annotated_calls <- left_join(non_annotated_calls, extra_annotations, by = c("patient_id" = "Sample", "CHROM", "POS", "REF", "ALT", "VAF", "DEPTH" = "Depth"))


rescued_calls <- 
  right_join(non_annotated_calls, rescue_mutations, by = c("CHROM" = "Chrom", "GENE" = "Gene", "REF" = "Ref", "POS" = "Position", "ALT" = "Alt"))


non_annotated_calls <- non_annotated_calls %>% 
  filter(!(VAF < 0.02 |
           # DEPTH < 50  |
           # (Left_reads + Right_reads) < 10 |
           (LOCATION == "intronic" | FUNCTION == "synonymous_SNV") |
           (dbSNP == "PRESENT" & VAF >= 0.450 & VAF <= 0.550) |
           # (dbSNP == "PRESENT" & VAF > 0.95) |
           COSMIC.y == "NOT_CONFIRMED_SOMATIC")
         )

write_csv(non_annotated_calls, paste0(path, "/05172021/Annotation Filtered Calls " , format(Sys.time(),"%m%d%Y"), ".csv"))


full_calls <- bind_rows(non_annotated_calls, rescued_calls) %>% 
  filter(!is.na(patient_id))

write_csv(non_annotated_calls, paste0(path, "/05172021/Annotated Rescued and Filtered Calls " , format(Sys.time(),"%m%d%Y"), ".csv"))














