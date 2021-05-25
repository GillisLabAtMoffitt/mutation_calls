# import library
library(tidyverse)
library(datapasta)

#######################################################################################  I  ### Load data----
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CIBMTR")
trusight <-
  read.delim(paste0(path, "/05172021/trusight-myeloid-amplicon-track.bed"))

panel <- tibble::tribble(
  ~`Gene;Target region (exon);Gene;Target region (exon)`,
                                   "ABL1; 4-6; KDM6A; full",
                          "ASXL1; 12; KIT; 2, 8-11, 13, 17",
          "ASXL2; 12, 13 (exon 13 CDS only); KLHDC8A; full",
                         "BCOR; full; KMT2A; 1, 3, 5-8, 27",
                                 "BCORL1; full; KRAS; 2, 3",
                               "BRAF; 15; MGA*; 3, 4, 6, 8",
                              "BRCC3; full; MPL; 4, 10, 12",
                                       "CALR; 9; MYB; full",
                                     "CBL; 8, 9; MYC*; 1-3",
                                  "CBLB; 9, 10; MYD88; 3-5",
                           "CBLC; 9, 10; NOTCH1; 26-28, 34",
                                   "CCND1*; 4, 5; NPM1; 11",
                                 "CCND2*; 4, 5; NRAS; 2, 3",
                         "CDKN2A; full; PDGFRA; 12, 14, 18",
                                 "CSF3R; 14-17; PHF6; full",
                "CUX1; full; PPM1D; 5, 6 (exon 6 CDS only)",
                                  "DHX15*; 3-5; PTEN; 5, 7",
                              "DNMT3A; ful;l PTPN11; 3, 13",
                              "ETV6/TEL; full; RAD21; full",
                                  "EZH2; full; RUNX1; full",
                         "FBXW7; 9-11; SETBP1; 4 (partial)",
                           "FLT3; 14, 15, 20; SF3B1; 13-16",
                           "GATA1; 2; SMC1A; 2, 11, 16, 17",
                 "GATA2; 2-6; SMC3; 10, 13, 19, 23, 25, 28",
                                    "GIGYF2*; 22; SRSF2; 1",
                                  "GNAS; 8, 9; STAG2; full",
                          "GNB1; 2, 3, 4, 5; STAT3; 20, 21",
                                   "HRAS; 2, 3; TET2; 3-11",
                                      "IDH1; 4; TP53; 2-11",
                                     "IDH2; 4; U2AF1; 2, 6",
                                   "IKZF1; full; WT1; 7, 9",
                           "JAK2; 12, 13, 14; ZBTB7A; full",
                                    "JAK3; 13; ZRSR2; full"
  ) %>% 
  separate(col = , into =  as.character(seq(1:4)), sep = "; ", extra = "warn", fill = "right")

a <- panel[c(1,2)] %>% `colnames<-`(c("Gene", "Target region (exon)"))
b <- panel[c(3,4)] %>% `colnames<-`(c("Gene", "Target region (exon)"))

panel <- rbind(a,b)
write_rds(panel, "TruSightMyeloid panel.rds")

panel <- readRDS("TruSightMyeloid panel.rds")
