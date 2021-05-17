# import library
library(tidyverse)


#######################################################################################  I  ### Load data----
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CIBMTR")

non_annotated_calls <-
  read_csv(paste0(path, "/05172021/CIBMTR_No Annotation Filtering_cleaned vcf.csv"))
extra_annotations <- 
  read.delim(paste0(path, "/05172021/MoffitCalls.annotated.txt"))

# clean extra to create var for join
extra_annotations <- extra_annotations %>% 
  mutate(CHROM = str_match(Variant, "chr(.*?):")[,2]) %>% 
  mutate(POS = as.numeric(str_match(Variant, ":([[:digit:]]*)")[,2])) %>% 
  mutate(REF = str_match(Variant, ":[[:digit:]]*(.*)>")[,2]) %>% 
  mutate(ALT = str_match(Variant, ">(.*)")[,2])

# Join with non_annotated_calls

non_annotated_calls <- left_join(non_annotated_calls, extra_annotations, by = c("patient_id" = "Sample", "CHROM", "POS", "REF", "ALT", "VAF", "DEPTH" = "Depth"))
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
