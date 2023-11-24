rm(list=ls())

options(scipen=999, # evitar notacion cientifica
        stringsAsFactors = FALSE, digits = 2)


# Cargar librerias
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(stargazer))
suppressMessages(library(corrr))
suppressMessages(library(zoo))
suppressMessages(library(haven))
suppressMessages(library(estimatr))
suppressMessages(library(lfe))
suppressMessages(library(haven))
suppressMessages(library(tidyverse))
suppressMessages(library(psych))
suppressMessages(library(ggcorrplot))
suppressMessages(library(sf))
suppressMessages(library(lavaan))

setwd("/Users/pedrofierroz/Dropbox/Academia/Projects/Paper_CJRES")

Data_2 <- read.csv('00_BBDD_2019_2021.csv', sep =";", header = T)


# SPATIAL INTERSECTION ----------------------------------------------------

Data_2_sp <- Data_2 %>% 
  filter(!is.na(Longitude))

Data_2_sp <- st_as_sf(Data_2_sp, coords = c('Longitude','Latitude'), crs = 4326, agr = 'identity')

ISMT <- st_read("/Users/pedrofierroz/Dropbox/Academia/Projects/Paper_CJRES/ISMT_Valpo.shp")

ISMT_GSE <- ISMT %>% 
  select(c(Alto, Medio, Bajo))

ISMT_Qi <- ISMT %>% 
  select(c(Q1, Q2, Q3, Q4, Q5))

ISMT_GSE$GSE <- colnames(ISMT_GSE)[apply(ISMT_GSE,1,which.max)]
ISMT_GSE$Qi <- colnames(ISMT_Qi)[apply(ISMT_Qi,1,which.max)]

summary(as.factor(ISMT_GSE$GSE))

Data_2_sp <- Data_2_sp %>%   #corrige vértices duplicados o geometrías "invalidas"
  st_join( # cruce de datos por operacion geografica
    ISMT_GSE,
    join = st_intersects, # asignaremos por interseccion
    left = FALSE # filtra valores donde no hubo interseccion
  ) 


Data_2_sp <- Data_2_sp %>% 
  mutate(
    Alto = ifelse(GSE=="Alto",1,0),
    Medio = ifelse(GSE=="Medio",1,0),
    Bajo = ifelse(GSE=="Bajo",1,0)
  )


###### TEST ########
model_B0 <- 'efin =~ efin2 + efin4 + efin5
          efin ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + periphery + camp'
fit_B0 <- sem(model = model_B0,
              #standardized = TRUE,
              data = Data_2)
pars.factors_B0 <- standardizedSolution(fit_B0)[ standardizedSolution(fit_B0)[,'op']=='~', c(3:7)]
pars.factors_B0

###### TEST ########
model_B00 <- 'efin =~ efin2 + efin4 + efin5
          efin ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + periphery'
fit_B00 <- sem(model = model_B00,
              #standardized = TRUE,
              data = Data_2)
pars.factors_B00 <- standardizedSolution(fit_B00)[ standardizedSolution(fit_B00)[,'op']=='~', c(3:7)]
pars.factors_B00



#MODELO CON GSE BAJO (A1 y A2)
model_A1 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + as.factor(Bajo) + periphery'

model_A2 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + as.factor(Bajo) + periphery + camp'

fit_A1 <- sem(model = model_A1,
              #standardized = TRUE,
              data = Data_2_sp)

fit_A2 <- sem(model = model_A2,
              #standardized = TRUE,
              data = Data_2_sp)

pars.factors_A1 <- standardizedSolution(fit_A1)[ standardizedSolution(fit_A1)[,'op']=='~', c(3:7)]
pars.factors_A2 <- standardizedSolution(fit_A2)[ standardizedSolution(fit_A2)[,'op']=='~', c(3:7)]

pars.factors_A1
pars.factors_A2

#MODELO CON I-MORAN p=0.01 (B1 y B2)
model_B1 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + periphery'

model_B2 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + periphery + camp'

fit_B1 <- sem(model = model_B1,
              #standardized = TRUE,
              data = Data_2)

fit_B2 <- sem(model = model_B2,
              #standardized = TRUE,
              data = Data_2)

pars.factors_B1 <- standardizedSolution(fit_B1)[ standardizedSolution(fit_B1)[,'op']=='~', c(3:7)]
pars.factors_B2 <- standardizedSolution(fit_B2)[ standardizedSolution(fit_B2)[,'op']=='~', c(3:7)]

pars.factors_B1
pars.factors_B2


## To extract RMSEA,SRMR and CFI
summary(fit_A1, fit.measures = TRUE)
summary(fit_A2, fit.measures = TRUE)


summary(fit_B1, fit.measures = TRUE)
summary(fit_B2, fit.measures = TRUE)


# To extract CFI
CFI_B1 <- fitMeasures(fit_B1, "cfi")
CFI_B1

CFI_B2 <- fitMeasures(fit_B2, "cfi")
CFI_B2


# CORR MATRIX -------------------------------------------------------------

Cor_m <- as.data.frame(Data_2_sp) %>% 
  select(GSE, Cluster_01,camp,educ,SES) %>% 
  filter(!is.na(SES)) %>% 
  filter(educ != "11") %>% 
  mutate(
    GSE = recode(GSE,
                 "Alto" = 0,"Medio" = 0,"Bajo" = 1),
    educ = recode(educ,
                  "Basica incompleta o menos"= 0,
                  "Basica completa"= 1,
                  "Media incompleta"= 2,
                  "Media completa"= 4,
                  "Media tecnica completa"= 3,
                  "Superior tecnica incompleta"= 5,
                  "Superior tecnica completa"= 6,
                  "Universitaria inompleta"= 7,
                  "Universitaria completa"= 8,
                  "Postgrado"= 9))


cor( Cor_m, method="spearman" )

# FE ----------------------------------------------------------------------


#REEMPLAZAMOS PERIFERIA POR EFECTOS DE COMUNA (multilcolineales)
model_B3 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + com1 + com2 + com3 + com4 + com5 + com6 + com7 + com8 + com9'

model_B4 <- 'efex =~ efex1 + efex3 + efex4
          efex ~ edad + SES + sexo + educ + yearcat3 + yearcat2 + Cluster_01 + camp + com1 + com2 + com3 + com4 + com5 + com6 + com7 + com8 + com9'

fit_B3 <- sem(model = model_B3,
              #standardized = TRUE,
              data = Data_2)

fit_B4 <- sem(model = model_B4,
              #standardized = TRUE,
              data = Data_2)

pars.factors_B3 <- standardizedSolution(fit_B3)[ standardizedSolution(fit_B3)[,'op']=='~', c(3:7)]
pars.factors_B4 <- standardizedSolution(fit_B4)[ standardizedSolution(fit_B4)[,'op']=='~', c(3:7)]

pars.factors_B3
pars.factors_B4