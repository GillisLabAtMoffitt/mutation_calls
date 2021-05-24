# import library
library(tidyverse)


#######################################################################################  I  ### Load data----
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CIBMTR")

extra_annotations <- 
  read.delim(paste0(path, "/05172021/MoffitCalls.annotated.txt"))

#######################################################################################  II  ### Clean data----

# 1.Clean extra_annotations to create var for join in step 5
extra_annotations <- extra_annotations %>% 
  mutate(CHROM = str_match(Variant, "chr(.*?):")[,2]) %>% 
  mutate(POS = as.numeric(str_match(Variant, ":([[:digit:]]*)")[,2])) %>% 
  mutate(REF = str_match(Variant, ":[[:digit:]]*(.*)>")[,2]) %>% 
  mutate(ALT = str_match(Variant, ">(.*)")[,2])

# write_csv(extra_annotations, paste0(path, "/05172021/MoffittCalls annotated separated " , format(Sys.time(),"%m%d%Y"), ".csv"))
# 2.Submit this file to Franklin

# 3.Prep rescued mutations from MoffitCalls.annotated analyzed by Franklin
franklin_mutations <-
  read.delim(paste0(path, "/05172021/export_1621355280244.tsv")) %>%
  # read.delim("/Users/colinccm/Documents/GitHub/Gillis/no vpn data/CIBMTR/05172021/export_1621355280244.tsv") %>%
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
  )), HgvsP1 = str_remove(HgvsP1, "p.")) %>%
#   # select("Chrom", "Position", "Ref", "Alt", "Gene", "HgvsC", "HgvsP", "Effect", "Region", "Genoox.Classification", "HgvsP1") %>% 
  mutate(Chrom = str_remove(Chrom, "chr"))

# 5.Join with non_annotated_calls with extra_annotations to get cosmic, right/left reads,...
extra_annotations <- right_join(extra_annotations,
                                 franklin_mutations,
                                 by = c("CHROM" = "Chrom", "POS" = "Position", 
                                        "REF" = "Ref", "ALT" = "Alt"))

# 4.Store the mutation in Franklin from our initial file which could be filtered out in the step 6
rescued_calls <- extra_annotations %>% 
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
  )

# 6.Filtering somatic mutations
extra_annotations <- extra_annotations %>% 
  filter(!(VAF < 0.02 |
           Depth < 50  | # NA
           (Left_reads + Right_reads) < 10 | # NA
           (Region == "INTRONIC" | Effect == "SYNONYMOUS") |
           (dbSNP == "PRESENT" & VAF >= 0.450 & VAF <= 0.550) |
           (dbSNP == "PRESENT" & VAF > 0.95) | # NA
           COSMIC == "NOT_CONFIRMED_SOMATIC")
         ) %>% 
  bind_rows(., rescued_calls)

write_csv(extra_annotations, paste0(path, "/05172021/Annotated Raw and Filtered Calls " , format(Sys.time(),"%m%d%Y"), ".csv"))














