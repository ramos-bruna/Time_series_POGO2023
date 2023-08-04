#Load package needed
library(lubridate)
library(plyr)
library(dplyr)
library(data.table)
library(imputeTS)
library(TTR)
library(vegan)
library(ade4)
library(factoextra)
library(tidyr)
library(reshape2)
library(gridExtra) #arrange plot
library(ggplot2) #plot
library(rshift)
library(funrar) #relative abundance

data(doubs)
?doubs
label = row.names(t(doubs$fish))
label
species = data.frame(doubs$fish)
env = data.frame(doubs$env)
plot(env)


#data exploration fish
head(species)
str(species)

#remove rare sp
#1st calculate total abundance
colnames(species)
mel_sp = reshape2::melt(species, measure.vars = c("Cogo","Satr", "Phph", "Neba", "Thth", "Teso",
                                                  "Chna", "Chto", "Lele", "Lece", "Baba", "Spbi", "Gogo",
                                                  "Eslu", "Pefl", "Rham", "Legi", "Scer", "Cyca", "Titi",
                                                  "Abbr", "Icme", "Acce", "Ruru", "Blbj", "Alal",
                                                  "Anan"), 
                        na.rm =T, value.name = "Abundance", variable.name = "Species")

mel_sp
tot_abund <-
  mel_sp %>%
  group_by(Species) %>%
  summarise(
    sum_abn = sum(Abundance, na.rm = F)
  )

print(tot_abund)

#2nd relative
# Calculate the total abundance for each sample (row) in the matrix
column_sum <- sum(tot_abund$sum_abn)

rel_abund <- (tot_abund$sum_abn / column_sum) *100
print(rel_abund)

#bind the name of the sp in the rel abund
sp = as.character(tot_abund$Species)
rel_abund_sp <- cbind(sp,rel_abund)

rel_abund_sp = as.data.frame(rel_abund_sp)
rel_abund_sp$V2 = as.numeric(rel_abund_sp$V2)
print(rel_abund_sp)

#looking all sp
ggplot(rel_abund_sp, aes(x = sp, y = V2)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Fishes Relative Abundance", x = "Species", y = "Relative Abundance (%)") +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) 


#Remove rare species
#We dont have rare species


#preliminary analysis
bartlett.test(species)

#Test the normality of each species:
for(i in 1:ncol(species)){
  qqnorm(species[,i],main = colnames(species[i]))
  qqline(species[,i], col = 2,)
}

#data transformation - in this case abundance
fish_std = decostand(species, "hellinger")
head(fish_std)

#fish PCA, we need for RDA after
PCA_fish = dudi.pca(fish_std, scale=T, center = T,nf=5) 
PCA_fish

s.corcircle(PCA_fish$co)

PCA_fish$co
PCA_fish$li

#Scree Plot, we see better the distribution in the dimentions
fviz_eig(PCA_fish, addlabels = TRUE)


# Graph of the species
fviz_pca_var(PCA_fish, col.var = "black")

# Result for species
res.fish = get_pca_var(PCA_fish)
res.fish
res.species$coord          # Coordinates
res.fish$contrib        # Contributions to axes

#SPECIES CONTRIBUTION
inertia.dudi(PCA_fish,col.inertia=T, row.inertia = T)


#We set a threshold of significant contribution 
# --> 1/(number of species - 1) *100
1/(nrow(PCA_fish$co))*100 #3.703704


#se tiver tempo aqui vale um subset com os dados >3.70
dados %>%
  filter(UR > 90)

#############################################################
#################starting to deal with env data
############################################################
install.packages("GGally")
library(GGally)

ggpairs(env)

env_data_long <- env %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

env_data_long

# Create a single plot with multiple facets (bar plot)
ggplot(env_data_long, aes(x = parameter, y = value)) +
  geom_bar(stat = "identity") +
  labs(title = "Environmental Parameters", x = "Parameter", y = "Value") +
  theme_minimal()


ggplot(env_data_long, aes(x = parameter, y = value)) +
  geom_boxplot() +
  labs(title = "Environmental Parameters", x = "Parameter", y = "Value") +
  theme_minimal()

ggplot(env_data_long, aes(x = NULL, y = value)) +
  geom_boxplot() +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(title = "Environmental Parameters", x = NULL, y = "Value") +
  theme_minimal()


#CENTERING AND SCALING
env_scale = as.data.frame(scale(env, center=TRUE, scale=TRUE))
nrow(env_scale)
nrow(fish_std)

# RDA with all the environmental variables (package vegan)
Fish_RDA = rda(fish_std,env_scale)
Fish_RDA #like this it env explain 70.81%

#test the signifnicace of the RDA
anova.cca(HR_RDA) #with vegan 
#it is significant

#test the colinearity between explanatory variables
#All number must be smaller than 10, otherwise remove one of the parameter
testVIF = vif.cca(Fish_RDA) 
testVIF

#Correlation between variables
lilibrary(corrplot)
Env_cor = cor(env, use = "pairwise.complete.obs")  #take the NA
Env_cor

Corrplot_env = corrplot(Env_cor, method = "number", type = "lower")
Corrplot_env

corrplot(Env_cor, method = "number", type = "lower", order = "hclust", tl.col = "black")

#removing column based on testVIF
colnames(env)
env_new = env[,-c(1,2,9,11)]
colnames(env_new)

#now everything again without dfs
env_new_CS = as.data.frame(scale(env_new, center=TRUE, scale=TRUE))

Fish_RDA_2 = rda(fish_std,env_new_CS)
Fish_RDA_2

#test the signifnicace of the RDA
anova.cca(Fish_RDA_2) #with vegan 
#it is significant

#test the colinearity between explanatory variables
#All number must be smaller than 10, otherwise remove one of the parameter
testVIF2 = vif.cca(Fish_RDA_2) 
testVIF2

#now everything again without dfs and amm 
colnames(env_new1)
env_new2 = env_new1[,-c(1)]

env_new_CS2 = as.data.frame(scale(env_new2, center=TRUE, scale=TRUE))

Fish_RDA_4 = rda(fish_std,env_new_CS2)
Fish_RDA_4

#test the signifnicace of the RDA
anova.cca(Fish_RDA_4) #with vegan 
#it is significant

#test the colinearity between explanatory variables
#All number must be smaller than 10, otherwise remove one of the parameter
testVIF4 = vif.cca(Fish_RDA_4) 
testVIF4

#now everything again without dfs,amm, alt and dbo
colnames(env_new2)
env_new3 = env_new2[,-c(8)]

env_new_CS3 = as.data.frame(scale(env_new3, center=TRUE, scale=TRUE))

Fish_RDA_4 = rda(fish_std,env_new_CS3)
Fish_RDA_4

#test the signifnicace of the RDA
anova.cca(Fish_RDA_4) #with vegan 
#it is significant

#test the colinearity between explanatory variables
#All number must be smaller than 10, otherwise remove one of the parameter
testVIF4 = vif.cca(Fish_RDA_4) 
testVIF4

############

#We applied a forward selection procedure to select a subset of 
#environmental variables to be used in the RDA
modcplt1 = rda(fish_std~.,env_new)

#Forward selection using ordistep() pack vegan
mod_empty = rda(fish_std~1, env_new3)

#variable selection with double stopping criteria R2 and pvalue
step.forward2 = ordiR2step(mod_empty, scope = formula(modcplt1),
                           direction = "both",
                           permutations = 999)

#Select the variables for the final:
select_var = env_scale[,c("oxy","flo","pho","nit","har")]


env_selected_fish = as.data.frame(scale(select_var, center=TRUE, scale=TRUE))

RDA_final = rda(fish_std,env_selected_fish)
RDA_final

anova.cca(RDA_final)

vif.cca(RDA_final)

#Adjusted R square
vegan::RsquareAdj(RDA_final)

summary(RDA_final)
ordiplot(RDA_final)


#############
#using ade4
#PCA
PCA_fish2 = dudi.pca(fish_std)

#RDA
RDA_ade4 = pcaiv(PCA_fish2,env_selected_fish)
summary(RDA_ade4)

plot(RDA_ade4)

#calculate with ADE4 the percentage of species explained by the environment:
sum(RDA_ade4$eig)/sum(PCA_fish2$eig)

#test the analysis. With ADE4 it is not possible to do an Anova like procedure.
#Instead we can do a permuation test
randtest(RDA_ade4)
plot(randtest(RDA_ade4))

#ask for only inertia of each canonical axis: contribution aux axes
inertia.dudi(RDA_ade4)
RDA_ade4$co
RDA_ade4$li

#Column (Species) absolute contributions
inertia.dudi(RDA_ade4, col = T, row = F)

(1/27)*100 #3.703704


RDA_ade4$cor #Correlation
RDA_ade4$fa #Canonical factor

#flo and nit influencing legi
2.567998e+00/2.336732e+00


####################
#plot RDA
#PLOT RDA
s.arrow(RDA_ade4$cor, xlim=c(-1,1), boxes = FALSE)
s.label(RDA_ade4$c1[-1,], add.plot=T, boxes = F)


# Extract species scores and site scores from the RDA result
species_scores <- scores(RDA_ade4, display = "species")
site_scores <- scores(RDA_ade4, display = "sites")

# Create a data frame for the biplot
biplot_data <- data.frame(
  x = c(site_scores$Axis1 , species_scores$Comp1),
  y = c(site_scores$Axis2, species_scores$Comp2),
  label = c(row.names(site_scores), row.names(species_scores)),
  type = c(rep("site", nrow(site_scores)), rep("species", nrow(species_scores)))
)

# Create the RDA biplot using ggplot2
ggplot(biplot_data, aes(x = x, y = y, label = label, color = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_text(size = 3, hjust = 0.5, vjust = 0.5) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  ggtitle("RDA Biplot")



plot(RDA_final, scaling = 1, display = c("cn","sp"), type = "points", pch = 20)
par(new=T)
plot(RDA_final, scaling = 1, display = c("cn","sp"), type = "text")


######



##########################################
#FIMMM

data(doubs)
pca1 <- dudi.pca(doubs$env, scan = FALSE)
pca2 <- dudi.pca(doubs$fish, scale = FALSE, scan = FALSE)
coiner1 <- coinertia(pca1, pca2, scan = FALSE)

if(adegraphicsLoaded()) {
  g1 <- s.corcircle(coiner1$aX, plot = FALSE)
  g2 <- s.value(doubs$xy, coiner1$lX[, 1], plot = FALSE)
  g3 <- s.value(doubs$xy, coiner1$lX[, 2], plot = FALSE)
  g4 <- s.arrow(coiner1$c1, plot = FALSE)
  g5 <- s.match(coiner1$mX, coiner1$mY, plot = FALSE)
  g6 <- s.corcircle(coiner1$aY, plot = FALSE)
  g7 <- s.arrow(coiner1$l1, plot = FALSE)
  g8 <- s.value(doubs$xy, coiner1$lY[, 1], plot = FALSE)
  g9 <- s.value(doubs$xy, coiner1$lY[, 2], plot = FALSE)
  G <- ADEgS(list(g1, g2, g3, g4, g5, g6, g7, g8, g9), layout = c(3, 3))
  
} else {  
  par(mfrow = c(3, 3))
  s.corcircle(coiner1$aX)
  s.value(doubs$xy, coiner1$lX[, 1])
  s.value(doubs$xy, coiner1$lX[, 2])
  s.arrow(coiner1$c1)
  s.match(coiner1$mX, coiner1$mY)
  s.corcircle(coiner1$aY)
  s.arrow(coiner1$l1)
  s.value(doubs$xy, coiner1$lY[, 1])
  s.value(doubs$xy, coiner1$lY[, 2])
  par(mfrow = c(1, 1))
}

adegraphicsLoaded
