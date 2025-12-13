

###### clean data weed seed
setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/plant_physiology")
weed <- read_csv("raw_weed_germination.csv")
head(weed)


f<-weed %>% filter(Species=="foxtail")
# tidy data  

colnames(f) <- c("Trt_ID" ,
                 "Treatment",
                 "Rep",
                 "Block",
                 "Species",
                 "foxtail_num_nongerm",
                 "foxtail_num_germ",
                 "foxtail_total_seeds",
                 "foxtail_prop_nongerm")
f$Species<-NULL
f$Block <-NULL
f$Treatment<-NULL
f$Rep<-NULL
head(f)
## pigweed


p<-weed %>% filter(Species=="pigweed")
# tidy data  

colnames(p) <- c("Trt_ID" ,
                 "Treatment",
                 "Rep",
                 "Block",
                 "Species",
                 "pigweed_num_nongerm",
                 "pigweed_num_germ",
                 "pigweed_total_seeds",
                 "pigweed_prop_nongerm")

p$Species<-NULL
weed<-full_join(p,f)
head(weed)
write.csv(weed, "weed_seed_decay.csv")
weed<-read.csv("weed_seed_decay.csv", row.names = 1)
head(weed)


