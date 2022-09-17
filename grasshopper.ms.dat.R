
# load packages
library(ggplot2)
library(car)
library(nlme)

# load data
g = read.csv("grasshopper.dat.5_13_2021.csv")

# clean data: seperate Melanoplus scudderi and remove intermediate ("low") glades
ms = g[which(g$species == "MS" & g$crot.pop != "low"),]


# linear mixed models  #############################################################################################

## dry mass ########################################################################################################

mod = lme(dry.mass ~ body.length + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(dry.mass ~ body.length + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

### explore and plot dry mass 
ggplot(data = ms, aes(y = dry.mass, x = body.length, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = dry.mass, x = body.length, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = dry.mass, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

## C:N ratio #######################################################################################################
mod = lme(C.N ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")  ## no sex:crot.pop interaction - remove in following model

mod2 = lme(C.N ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

### explore and plot C:N
ggplot(data = ms, aes(y = C.N, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = C.N, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = C.N, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

### C percent ######################################################################################################

mod = lme(C.percent ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(C.percent ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

ggplot(data = ms ,aes(y = C.percent, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms,aes(y = C.percent, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = C.percent, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

### N percent ######################################################################################################
mod = lme(N.percent ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(N.percent ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")
lsmeans(mod2, pairwise ~  sex, adjust="tukey") # post hoc analysis of sex effect

ggplot(data = ms, aes(y = N.percent, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = N.percent, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = N.percent, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

## lipid:protein  ##################################################################################################

mod = lme(lipid.protein ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod = lme(lipid.protein ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

ggplot(data = ms, aes(y = lipid.protein, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = lipid.protein, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = lipid.protein, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()


### lipid percent ##################################################################################################

ms2 = ms[which(ms$lipid.percent > 0),] # remove negative lipid.percent values 

mod = lme(lipid.percent ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms2, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(lipid.percent ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms2, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")
lsmeans(mod2, pairwise ~ crot.pop, adjust="tukey")

ggplot(data = ms2, aes(y = lipid.percent, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms2, aes(y = lipid.percent, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms2, aes(y = lipid.percent, x = sex))+
  geom_boxplot()+
  theme_classic()

ggplot(data = ms, aes(y = lipid.percent, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()


### protein percent ################################################################################################

mod = lme(protein.percent ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(protein.percent ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")
lsmeans(mod2, pairwise ~ crot.pop, adjust="tukey") 

ggplot(data = ms, aes(y = protein.percent, x = sex))+
  geom_boxplot()+
  theme_classic()

ggplot(data = ms, aes(y = protein.percent, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = protein.percent, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = protein.percent, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()


### C13 ###########################################################################################################
mod = lme(C13 ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(C13 ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

ggplot(data = ms, aes(y = C13, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = C13, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = C13, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

### N15 ###########################################################################################################
mod = lme(N15 ~ dry.mass + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ts, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(N15 ~ dry.mass + sex + crot.pop, random = ~1|glade, data = ts, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

ggplot(data = ms, aes(y = N15, x = dry.mass, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = N15, x = dry.mass, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = N15, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()

### Body Length ########################################################################################################################################
mod = lme(body.length ~ sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod2 = lme(body.length ~ sex + crot.pop, random = ~1|glade, data = ms, method = "REML")
anova.lme(mod2, type = "marginal")
Anova(mod2, type = "III")

ggplot(data = ms, aes(y = body.length, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()


# body morphology PCA  ###########################################################################################

# prepare new data table
tibia = rowMeans(cbind(ms$tibia.length.l, ms$tibia.length.r), na.rm = T) # average left and right leg measurements
femur = rowMeans(cbind(ms$femur.length.l, ms$femur.length.r), na.rm = T)
femur.w = rowMeans(cbind(ms$femur.width.l, ms$femur.width.r), na.rm = T)

ms.body = cbind(ms[,c("glade","crot.pop","sex","body.length")], tibia, femur, femur.w)
ms.body = na.omit(ms.body) # remove NAs for prcomp()

# check assumption of correlation for PCA
pairs(ms.body[,c("body.length","tibia","femur","femur.w")])
cor(ms.body[,c("body.length","tibia","femur","femur.w")])


# perform pca 
morph.pca = prcomp(ms.body[,c("body.length","tibia","femur","femur.w")])
summary(morph.pca) # PC1 accounts for most of the variation, only factor used in LMMs
morph.pca$x[1:2] 
morph.pca$sdev^2
morph.pca$rotation

ms.pca = cbind(ms.body, morph.pca$x) # add PCA values to ms.body

# linear mixed models
mod = lme(PC1 ~ body.length + sex + crot.pop + sex:crot.pop, random = ~1|glade, data = ms.pca, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

mod = lme(PC1 ~ body.length + sex + crot.pop, random = ~1|glade, data = ms.pca, method = "REML")
anova.lme(mod, type = "marginal")
Anova(mod, type = "III")

ggplot(data = ms.pca, aes(y = PC1, x = body.length, col = crot.pop))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms.pca, aes(y = PC1, x = body.length, col = sex))+
  geom_smooth(method = "lm",se = T)+
  geom_point()+
  theme_classic()

ggplot(data = ms, aes(y = PC1, x = sex, col = crot.pop))+
  geom_boxplot()+
  theme_classic()
