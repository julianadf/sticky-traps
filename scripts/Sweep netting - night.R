rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(lattice)
library(ggpubr)
library(effects)
library(visreg)
library(DHARMa)
library(MASS)
library(lme4)
library(car)
library(png)
library(cowplot)
library(visreg)
library(pscl)
library(glmmTMB)

# Load and check data ----
sticky.traps.raw <- read.csv2(here("data", "Sticky_traps_raw.csv" ))
day.counts.raw <- read.csv2(here("data", "Abundance day-flying pollinators.csv" ))
night.counts.raw <- read.csv2(here("data", "SweepNets.csv" ))
site.data <- read.csv2(here("data", "Site metadata.csv" ))

# Night-flying pollinator analysis ----

# Check and transform night insect data (ALL INSECTS - spiders)
night.counts <- night.counts.raw %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(Mown.side = as.factor(Mown.side)) %>% 
  mutate(Mown.site = as.factor(Mown.site)) %>% 
  mutate(Round = as.factor(Round)) %>% 
  #mutate(Others=as.integer(Others)) %>% 
  group_by(Site, Mown.site, Round) %>% 
  summarise_at(vars(Staphylinidae, Cantharidae, Coccinellidae, Curculionidae, pollen.beetles, flea.beetles,
                    other.Coleoptera, Parasitic.wasps, Formicidae, Other.Apocrita, Symphyta, Diptera..Syrphidae, 
                    Diptera..other.Brachycera, Diptera.other.Nematocera, Lepidoptera..butterfly.adults, 
                    Lepidoptera..butterfly.larvae, Lepidoptera..moth.adults, Lepidoptera..moth.larvae, 
                    Aphids, Cicadellidae, Delphacidae, Other.Homoptera, Pentatomidae, Miridae, other.Heteroptera, 
                    Neuroptera, Thysanoptera, Collembola, Dermpatera, Orthoptera, Acari), .funs=sum)
night.counts

# Join day-flying data with site info & turn it into long format
nightverge.database <-left_join(night.counts, site.data, by="Site")
night.flying <- gather(nightverge.database, insect.group, indiv, Staphylinidae:Acari, factor_key = TRUE)

ggplot(night.flying, aes(x=log10(indiv))) + geom_histogram(fill="lightgreen", color="grey50")

# model attempt 1: 
test.model <- glmer(indiv~ Road.verge * scale(Traffic) + Round + (1|insect.group), offset=log10(Verge.width),
                     family="poisson", data = night.flying)
summary(test.model)
plot(allEffects(test.model))
mod_dharma1 <- test.model %>% simulateResiduals(n=1000)
plot(mod_dharma1)
mod_dharma1 %>% testZeroInflation() # data is zero inflated
mod_dharma1 %>% testDispersion()

# model attempt 2: zero-inflated GLMM poisson
test.model2 <- formula(indiv~ Road.verge * scale(Traffic) + Round + (1|insect.group), 
                       offset=log10(Verge.width))
zip1 <- zeroinfl(test.model2, dist="poisson", link="logit", data=night.flying)
