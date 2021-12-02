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

# Load and check data ----
sticky.traps.raw <- read.csv2(here("data", "Sticky_traps_raw.csv" ))
day.counts.raw <- read.csv2(here("data", "Abundance day-flying pollinators.csv" ))
night.counts.raw <- read.csv2(here("data", "SweepNets.csv" ))
site.data <- read.csv2(here("data", "Site metadata.csv" ))

# Day-flying pollinator analysis ----

# Check and transform Day-flying insect data
day.counts <- day.counts.raw %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(Mown.side = as.factor(Mown.side)) %>% 
  mutate(Mown.site = as.factor(Mown.site)) %>% 
  mutate(Round = as.factor(Round)) %>% 
  mutate(Observer = as.factor(Observer)) %>% 
  group_by(Site, Mown.site, Round, Observer) %>% 
  summarise_at(vars(Butterflies, Bumblebees, Hoverflies, Bees), .funs=sum)
day.counts

# Join day-flying data with site info & turn it into long format
dayverge.database <-left_join(day.counts, site.data, by="Site")
day.flying <- gather(dayverge.database, insect.group, indiv, Butterflies:Bees, factor_key = TRUE)

ggplot(day.flying, aes(x=Traffic, y=indiv, colour=insect.group)) + geom_point()
hist(day.flying$indiv, breaks = 20)

#model attempt 1: don't have enough data
day.flying$obs <- c(1:240)
test.model <- glmer(indiv~ Road.verge * scale(Traffic) + Round + (1|insect.group),  offset=log10(Verge.width),
                    family="poisson", data = day.flying)
summary(test.model)
plot(allEffects(test.model))
mod_dharma1 <- test.model %>% simulateResiduals(n=1000)
plot(mod_dharma1)
mod_dharma1 %>% testZeroInflation()
mod_dharma1 %>% testDispersion()

# model attempt 2:
day.counts.aggr <- day.counts.raw %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(Mown.side = as.factor(Mown.side)) %>% 
  mutate(Mown.site = as.factor(Mown.site)) %>% 
  mutate(Round = as.factor(Round)) %>% 
  mutate(Observer = as.factor(Observer)) %>% 
  group_by(Site) %>% 
  summarise_at(vars(Butterflies, Bumblebees, Hoverflies, Bees), .funs=sum)
day.counts.aggr

# Join day-flying data with site info & turn it into long format
dayverge.database.aggr <-left_join(day.counts.aggr, site.data, by="Site")
day.flying.aggr <- gather(dayverge.database.aggr, insect.group, indiv, Butterflies:Bees, factor_key = TRUE)

ggplot(day.flying.aggr, aes(x=Traffic, y=indiv, colour=insect.group)) + geom_point()

day.flying.aggr$obs <- c(1:80)
test.model2 <- glmer(indiv~ Road.verge + scale(Traffic) + (1| insect.group),  family="poisson", data = day.flying.aggr)
summary(test.model2)
vif(test.model2)
plot(allEffects(test.model2))
#visreg(test.model2, "insect.group", by="Traffic", scale="response")
mod_dharma1 <- test.model2 %>% simulateResiduals(n=1000)
plot(mod_dharma1)

# model attempt 3:

# Check and transform Day-flying insect data
day.counts.allins <- day.counts %>% 
  mutate(tot.abund=sum(Butterflies, Bumblebees, Hoverflies, Bees))
day.counts.allins

# Join day-flying data with site info
day.counts.allins <-left_join(day.counts.allins, site.data, by="Site")

#model attempt 3: don't have enough data
day.counts.allins$obs <- c(1:60)
test.model3 <- glm(tot.abund ~ Road.verge + scale(Traffic) + Round, 
                     family="quasipoisson", data = day.counts.allins)
summary(test.model3)
plot(allEffects(test.model3))
mod_dharma1 <- test.model3 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
mod_dharma1 %>% testDispersion()

#model attempt 4: 
day.counts.allaggr <- day.counts.raw %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(Mown.side = as.factor(Mown.side)) %>% 
  mutate(Mown.site = as.factor(Mown.site)) %>% 
  mutate(Round = as.factor(Round)) %>% 
  mutate(Observer = as.factor(Observer)) %>% 
  group_by(Site) %>% 
  summarise_at(vars(Butterflies, Bumblebees, Hoverflies, Bees), .funs=sum)
day.counts.allaggr

test.model3 <- glm(tot.abund ~ Road.verge + scale(Traffic) + Round, 
                   family="quasipoisson", data = day.counts.allins)
summary(test.model3)
plot(allEffects(test.model3))
mod_dharma1 <- test.model3 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
mod_dharma1 %>% testDispersion()
mod_dharma1 %>% testZeroInflation()

# Final model?
test.model4 <- glmer.nb(indiv~ Road.verge * scale(Traffic) + Mown.site + (1|insect.group),  offset=log10(Verge.width),
                       data = day.flying)
summary(test.model4)
plot(allEffects(test.model4))
mod_dharma1 <- test.model4 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
mod_dharma1 %>% testDispersion()
mod_dharma1 %>% testZeroInflation()


# Check and transform sticky traps data
sticky.traps <- sticky.traps.raw %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(Road.verge=as.factor(Road.verge)) %>% 
  mutate(Round=as.factor(Round)) %>% 
  mutate(Day.Night=as.factor(Day.Night)) %>% 
  group_by(Site, Round, Day.Night, Road.verge) %>%    
  summarise_at(vars(-Traffic, -Verge.width, -Serial.number, -Date, -Time, -Others), .funs=sum) 
sticky.traps



