##### Haematopus palliatus Trophic Ecology #####
rm(list=ls())
piru.piru <- read.csv (file="siber_local_temporada_names.csv", h=T, sep=";")
##### Kruskal Wallis com posthoc de Dunn entre consumidores####
### para d13C ##
library(FSA)
kruskal.test(iso1~group)

# Post hoc test (dunnTest)
PHC <- dunnTest(iso1~group,
                data=piru.piru,
                method="bh")
PHC 

## para d15N##
kruskal.test(iso2~group)
# Post hoc test (dunnTest)
PHN <- dunnTest(iso2~group,
                data=dados,
                method="bh")
PHN
##### Kruskall Wallis entre os anos ou sexo em cada local ######
rm(list=ls())
dados <- read.csv(file = "siber_sexo_itapeva.csv", h=T, sep=";") 

attach(dados)
boxplot(iso1 ~ group) 
boxplot(iso2 ~ group) 

# Kruskal Wallis para d13C
kruskal.test(iso1~group)
library(FSA)

# Post hoc test (dunnTest)
PHC <- dunnTest(iso1~group,
                data= dados,
                method="bh")
PHC 

kruskal.test(iso2~group)
library(FSA)
# Post hoc test (dunnTest)
PHN <- dunnTest(iso2~group,
                data=dados,
                method="bh")
PHN

##### SIBER ####
install.packages("rjags") 
library(rjags)
library(jagsUI)
library(coda)
library(ggplot2)
library(sp)
library(splancs)
library(SIBER)
library(hrbrthemes)
library(ellipse)
library(dplyr)
#remove items and remove graphs
rm(list=ls())
graphics.off()

set.seed(1) #always at the begining

###### SIBER TODOS OS LOCAIS ####
library(SIBER)
library(dplyr)
library(ggplot2)

setwd("C:/Users/gtnbi/Desktop/lais/")
ita.sex<-read.csv(file="siber_sexo_itapeva.csv", sep=";", h=T)

data <- ita.sex %>% mutate(group = factor(group), 
                           community = factor(community),
                           d13C = iso1, 
                           d15N = iso2,
                           .keep = "unused") 

# SIBER plot
plot <- ggplot(data = data, 
               aes(x = d13C, 
                   y = d15N)) + 
  geom_point(aes(color = group), size = 2) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 12),
        text = element_text(size = 14)) +
  scale_color_viridis_d()


# Including ellipses of 0.95
p.ell <- 0.95

x_limits <- c(-17.5, -11.5)
y_limits <- c(9, 16)


ellipse.plot <- plot +
  stat_ellipse(aes(group = interaction(group, community),
                   fill = group,
                   color = group),
               alpha = 0.25,
               level = 0.95,
               type = "norm",
               geom = "polygon") +
  labs(title = "Itapeva")+
  scale_color_brewer(palette = "Paired") +
  scale_fill_viridis_d() +
  coord_equal(xlim = x_limits, ylim = y_limits)+
  theme(plot.title = element_text(hjust = 0.5))


ellipse.plot

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siberobject)
print(group.ML)
# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics
# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siberobject) 
print(community.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siberobject, parms, priors)

## ----create-ellipse-df for ggplot---------------------------------

# how many of the posterior draws do you want?
n.posts <- 5

# decide how big an ellipse you want to draw
p.ell <- 0.95

# for a standard ellipse use
p.ell <- pchisq(1,2)

# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")

# now we need the group and community names

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

first.plot <- ggplot(data = dados, aes(iso1, iso2)) +
  geom_point(aes(color = factor(group)), size = 1)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=15)) 
print(first.plot)

second.plot <- first.plot + facet_wrap(~factor(group))
print(second.plot)

third.plot <- second.plot + 
  geom_polygon(data = dados,
               mapping = aes(iso1, iso2,
                             group = group,
                             color = factor(group),
                             fill = NULL),
               fill = NA,
               alpha = 0.2)
print(third.plot)

plotSiberObject(siberobject,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
                cex=1.0,
                y.limits = c(9,16),
                x.limits = c(-17,-12),
                main= "")

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)
colnames(group.ML) <- c("PT", "PG", "ITA", "PC", "LP")

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML),
                 xlab="",
                 ylab = "",
                 bty = "L",
                 las = 2
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

##### SIBER ellipses overlap between sites ####

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "1.Itapeva" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "1.Passo de Torres"

# The overlap of the maximum likelihood fitted standard ellipses are 
# estimated using
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siberobject, 
                             p.interval = NULL, n = 100)

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siberobject, 
                                   p.interval = 0.95, n = 100)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])

# itapeva e passo de torres
ellipse1 <- "1.Itapeva" 
ellipse2 <- "1.Passo de Torres"
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siberobject, 
                             p.interval = NULL, n = 100)
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, siberobject, 
                                   p.interval = 0.95, n = 100)

prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] + 
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])


##### SIBER Lagoa do Peixe MACHOS E FEMEAS ######
dados <- read.csv(file="siber_sexo_lagoadopeixe.csv", sep=";" , header=T) # load data

siberobject <- createSiberObject(dados) # create siber object 

#create arguments separated from the original siber function #

community.hulls.args <- list(palette(c("magenta", "blue")), 
                             lty = 1, lwd = 1)

group.ellipses.args  <- list(n = 10000, p.interval = 0.95, lty = 1, lwd = 4)

plot(dados$iso1,dados$iso2,xlab=expression({delta}^13*C~'\u2030'),
     ylab=expression({delta}^15*N~'\u2030'),
     xlim=c(-15,-13),ylim=c(9,16),
     xaxt="n", yaxt="n", main="Lagoa do Peixe")
axis(1, at = seq(-15, -13, by = 1), las=2)
axis(2, at = seq(9, 16, by = 1), las=2)

points(dados$iso1,dados$iso2,pch=c(21,21,21)[dados$group], ##nao esta plotando os pontos
       bg=c("magenta","blue")[dados$group], cex= 1.2)

legend("bottomright", c("Machos", "Femeas"),
       fill=c("blue", "magenta"), lty=0, bty="n")

plotGroupEllipses(siberobject, n = 10000, p.interval = 0.95,
                  lty = 1, lwd = 2)


# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siberobject)
print(group.ML)

# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics

# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siberobject) 
print(community.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siberobject, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Group"), 
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses Lagoa do Peixe"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

##### SIBER Praia Grande MACHOS E FEMEAS ######
dados <- read.csv(file="siber_sexo_pg.csv", sep=";" , header=T)

siberobject <- createSiberObject(dados)

community.hulls.args <- list(palette(c("magenta", "blue")), 
                             lty = 1, lwd = 1)

group.ellipses.args  <- list(n = 10000, p.interval = 0.95, lty = 1, lwd = 4)

plot(dados$iso1,dados$iso2,xlab=expression({delta}^13*C~'\u2030'),
     ylab=expression({delta}^15*N~'\u2030'),
     xlim=c(-18,-11),ylim=c(7,16),
     xaxt="n", yaxt="n", main="Praia Grande")
axis(1, at = seq(-18, -11, by = 2), las=2)
axis(2, at = seq(7, 16, by = 2), las=2)

points(dados$iso1,dados$iso2,pch=c(21,21)[dados$group], 
       bg=c("magenta","blue")[dados$group], cex= 1.2)

legend("bottomright", c("Machos", "Femeas"),
       fill=c("blue", "magenta"), lty=0, bty="n")

plotGroupEllipses(siberobject, n = 10000, p.interval = 0.95,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siberobject)
print(group.ML)
# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics
# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siberobject) 
print(community.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siberobject, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)


siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Group"), 
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses Praia Grande"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

###### SIBER VARIACAO ANUAL NA PRAIA GRANDE ####
dados <- read.csv(file="siber_anos_pg.csv", sep=";" , header=T)

siberobject <- createSiberObject(dados)

community.hulls.args <- list(palette(c("magenta", "blue", "yellow", "orange")), 
                             lty = 1, lwd = 1)

group.ellipses.args  <- list(n = 10000, p.interval = 0.95, lty = 1, lwd = 4)
par(mfrow=c(1,1))
plot(dados$iso1,dados$iso2,xlab=expression({delta}^13*C~'\u2030'),
     ylab=expression({delta}^15*N~'\u2030'),
     xlim=c(-18,-11),ylim=c(7,16),
     xaxt="n", yaxt="n", main="Lagoa do Peixe")
axis(1, at = seq(-18, -11, by = 2), las=2)
axis(2, at = seq(7, 16, by = 2), las=2)

points(dados$iso1,dados$iso2,pch=c(21,21,21,21)[dados$group], 
       bg=c("magenta","blue", "yellow", "orange")[dados$group], cex= 1.2)

legend("bottomright", c("2018", "2019", "2020", "2021"),
       fill=c("blue", "magenta", "yellow", "orange"), lty=0, bty="n")

plotGroupEllipses(siberobject, n = 10000, p.interval = 0.95,
                  lty = 1, lwd = 2)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siberobject)
print(group.ML)
# A second plot provides information more suitable to comparing
# the two communities based on the community-level Layman metrics
# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siberobject) 
print(community.ML)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siberobject, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Group"), 
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)


#### medias e SD isotopos presas e consumidores para modelos de mistura####
rm(list=ls())
graphics.off()

dados <- read.csv(file="2023 isotopos_presas.csv", h=T, sep = ";")

mean_d15N <- aggregate(dados$d15N, by= list(dados$species), FUN=mean)
mean_d13C <- aggregate(dados$d13C, by= list(dados$species), FUN=mean)
sd_d15N <- aggregate(dados$d15N, by= list(dados$species), FUN=sd)
sd_d13C <- aggregate(dados$d13C, by= list(dados$species), FUN=sd)

presas <- cbind(mean_d13C, mean_d15N$x, sd_d13C$x, sd_d15N$x)
colnames(presas) <- c("Presas", "mean_d13C", "mean_d15N", "sd_d13C", "sd_d15N")
presas

## media e desvio dos dados de piru-piru
mean_d15N <- aggregate(dados$iso2, by= list(dados$group), FUN=mean)
mean_d13C <- aggregate(dados$iso1, by= list(dados$group), FUN=mean)
sd_d15N <- aggregate(dados$iso2, by= list(dados$group), FUN=sd)
sd_d13C <- aggregate(dados$iso1, by= list(dados$group), FUN=sd)

#### SIMMR MODELOS DE MISTURA####
# remove previously loaded items from the current environment and remove previous graphics.
rm(list=ls())
graphics.off()

library(vegan)
library(simmr)

##### SIMMR modelo com dados de 100 consumidores, locais PT-9, PG-42, ITA-18, CB- 7, LP-24 ####


mix <- matrix(c(-13.051768, -13.2236656, -13.0476752, -13.03642, -13.0118632, -13.2021784, -13.2594776, -14.1198085,
                -14.38885, -16.3829466, -13.8294, -14.7768832, -14.330768, -14.1199888, -15.8021972, -13.7256085,
                -13.7485672, -13.2635704, -13.435468, -13.5920176, -13.1899, -13.4620712, -13.742428, -14.0094832,
                -13.6646648, -13.819168, -13.3423568, -13.8590728, -13.686152, -15.7478545, -15.084613, -15.291568,
                -15.973534, -15.5822905, -16.1420545, -15.055048, -14.9180635, -15.210757, -14.4253135, -14.773195,
                -14.434183, -14.9963, -15.1633, -15.2433, -14.6367, -14.9222, -15.0408, -14.9341, -14.7898, -14.6278,
                -15.5437, -15.3016, -15.2058, -15.0151, -14.7552, -15.1959, -15.4043, -15.2423, -15.3046, -15.0931,
                -13.067116, -13.2072944, -13.6646648, -14.276503, -14.745601, -15.33493, -14.7544705, -14.869774, -15.068845,
                -14.0192875,  -14.319865, -14.2873435, -14.0468815, -14.7544705,  -13.9286215, -14.258764, -13.693087, -13.866535,
                -14.01436, -13.8497815, -13.9365055, -13.622131, -14.19175, -14.1730255, -14.0409685, -13.8103615, -14.4607915, -14.258764,
                -14.138533, -14.0961565, -13.894129, -13.8364, -13.6892, -13.7139, -13.6082, -13.7336, -13.5864, -13.4837, -13.8364, -14.1644,
                13.1703546, 12.5025446,12.6545998,13.1035736,13.0593954,12.8724086,13.5083692,14.9667616,14.2428658,9.2196328,12.750148,
                12.3515168,12.4748048,12.3936402,8.9641912,14.1121912, 14.3076864,14.1915902,13.7385068,13.3388482,13.5998078,13.8104248,
                13.4159032,12.755285,13.1847382,12.4295992, 13.063505, 11.799803,12.1306258, 12.1136995, 12.8002597, 12.3242308, 12.8396695,
                13.2119884,11.2622404, 11.0807479, 11.2197193, 12.2889694,14.0323345,14.5519216, 12.4331263, 13.58261, 13.19701,12.99654,13.49158,
                13.107, 12.99245, 12.972, 12.32865, 10.83741, 13.19394,13.97741,13.52635, 14.21777,13.99991, 13.47214,13.76671,13.70739, 13.58261,
                13.99275,14.034398, 14.096042, 13.8063152, 14.8682371, 14.0976718, 13.8072838, 14.1899737, 12.8770051, 12.1894078, 13.1279833, 12.7587757,
                13.0968703,11.6013721, 13.7243158, 13.8622501, 13.2949564, 11.7590113, 11.3317261, 11.7818275, 12.3273421, 11.9861362, 12.388531,
                12.7774435, 11.880352, 12.0856978, 12.0753268, 11.5868527, 11.4292135, 12.0307315,10.9438507,11.7268612, 13.00063, 13.3402, 13.07632, 13.107,
                12.75925, 13.27065, 12.97813, 12.82471, 12.89835), ncol=2, nrow=100)
colnames(mix) <-  c("d13C", "d15N")

## agrupando os locais 1-PT, 2-PG, 3-ITA, 4-Cabras, 5-PNLP
grp <- as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 3, 3, 3,
                    4, 4, 4, 4, 4, 4, 4,
                    5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                    5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                    5, 5, 5, 5))

#dados de todas as 11 presas para rodar biplot de carbono e nitrogenio #

s_names <- c("A_mactroides", "C_sapidus", "D_hanleyanus", "E_brasiliensis",
             "Isopodes", "N_granulata", "O_auriculata", "P_perna", "Polichaeta",
             "Stramonita", "Tagellus")
s_means <- matrix(c(-14.34499,-17.37590,-14.45716,-13.94067, -13.75714, -19.67770, -12.97285, -14.45223, -12.93064, -13.80152, -14.25975,    
                    10.53560, 11.26812, 10.49058, 10.44184,  11.41262,  10.32505,  13.87949,  10.22637,  10.47439,  11.58739,  10.12610), ncol=2, nrow=11) #primeiro carbono depois nitrogenio
s_sds <- matrix(c(0.85304723, 1.41502461, 0.64732883, 0.18102714, 0.00557483, 0.74084194,0.85554136, 0.73409650, 1.07772855, 0.48665711, 0.65922363,    
                  0.5478120,  1.0829590,  1.1503264,  0.7023311,  0.2082687,  1.6203752, 0.6918654,  1.4103151,  2.1370971, 0.8873974, 1.1388777), ncol=2, nrow=11)

# TEF de 11 presas para serem adicionados em cada presa#

c_means <- matrix(c(0.2,	0.2,	0.2, 0.2,	0.2, 0.2, 0.2, 0.2,	0.2, 0.2, 0.2,
                    2.7,	2.7,  2.7, 2.7,	2.7, 2.7, 2.7,  2.7, 2.7,	2.7, 2.7 ), ncol=2, nrow=11)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4,	0.4, 0.4, 0.4,	0.4, 0.4,
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,	0.4, 0.4), ncol=2, nrow=11)

##### SIMMR Codigo PRINCIPAL para carregar as matrizes dentro do simmr #####
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))

#GRAFICA OS DADOS: BIPLOT
plot(simmr_groups,group=1:5,xlab=expression(paste(delta^13, "C (\u2030)",sep="")),ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Diet of the American Oystercatcher',mix_name = NULL)
##### SIMMR para PASSO DE TORRES, PRAIA GRANDE e ITAPEVA ####
# rodar consumidores do primeiro codigo
# carregando dados das presas
s_names <- c("Amarillodesma", "Donax", "Emerita","Excirolana","Olivancillaria", "Perna", "Polichaeta")
s_means <- matrix(c(-14.34499,-14.45716,-13.94067,-13.75714, -12.97285, -14.45223, -12.93064,
                    10.53560, 10.49058, 10.44184, 11.41262, 13.87949, 10.22637, 10.47439), ncol=2, nrow=7) #primeiro carbono depois nitrogenio
s_sds <- matrix(c(0.85304723, 0.64732883,0.18102714, 0.00557483, 0.85554136, 0.73409650, 1.07772855,
                  0.5478120,1.1503264, 0.7023311, 0.2082687, 0.6918654, 1.4103151, 2.1370971 ), ncol=2, nrow=7)

#TEF
c_means <- matrix(c(0.2,	0.2,	0.2,0.2,	0.2,	0.2, 0.2, 
                    2.7,	2.7, 2.7,2.7,	2.7, 2.7, 2.7 ), ncol=2, nrow=7)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4,	0.4,
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=7)


#CARREGA AS MATRIZES NO SIMMR
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))

simmr_out <- simmr_mcmc(simmr_groups)

# Diagnostico de convergencia
summary(simmr_out,type='diagnostics')

#AVALIANDO O AJUSTE DO MODELO
posterior_predictive(simmr_out)

#Comparando as distribuicoes priori e posteriori do modelo Bayesiano
# Note que as curvas de probabilidade a priori - para prioris nao-informativas -
# consideram que todas as presas teriam iguais contribuicoes (1/5) para a dieta do consumidor
prior_viz(simmr_out, group=1)
prior_viz(simmr_out, group=2)
prior_viz(simmr_out, group=4)

# Media e desvios das posterioris
# Utilizei a funcao c() para mostrar os dois grupos ao mesmo
# tempo pois pedindo separadamente o R so estava dando as 
# estimativas do grupo 1
summary(simmr_out,type='statistics', group=c(1, 2,3,4))

# Quantis das posterioris. 
summary(simmr_out,type='quantiles', group = c(1, 2, 3,4))


compare_sources(simmr_out,source_names=c("Amarillodesma", "Donax", "Emerita","Excirolana","Olivancillaria", "Perna", "Polichaeta"), gr= 1)
compare_sources(simmr_out,source_names=c("Amarillodesma", "Donax", "Emerita","Excirolana","Olivancillaria", "Perna", "Polichaeta"), gr= 2)
compare_sources(simmr_out,source_names=c("Amarillodesma", "Donax", "Emerita","Excirolana","Olivancillaria", "Perna", "Polichaeta"), gr= 3)


plot(simmr_out,type='density', group = 1, title='Passo de Torres')
plot(simmr_out,type='density', group = 2, title='Praia Grande')
plot(simmr_out,type='density', group = 3, title='Itapeva')

plot(simmr_out,type='boxplot', group = 4, title='Praia Grande 2021')

plot(simmr_out,type='boxplot', group = 2, title='Praia Grande')

plot(simmr_out,type='boxplot', group = 3, title='Itapeva')

##### SIMMR para PRAIA DAS CABRAS #####
s_names <- c("Amarillodesma", "Donax", "Emerita", "Excirolana", "Polichaeta", "Ollivancillaria")
s_means <- matrix(c(-14.34499, -14.45716, -13.94067, -13.75714, -12.93064, -12.97285,   
                    10.53560,   10.49058, 10.44184,   11.41262, 10.47439, 13.87949 ), ncol=2, nrow=6) 
s_sds <- matrix(c(0.85304723, 0.64732883, 0.18102714, 0.00557483, 1.07772855, 0.85554136,
                  0.5478120,  1.1503264,  0.7023311,  0.2082687, 2.1370971, 0.6918654), ncol=2, nrow=6)

#TEF

c_means <- matrix(c(0.2,	0.2,	0.2, 0.2,	0.2, 0.2,
                    2.7,	  2.7,  2.7, 2.7,	2.7, 2.7), ncol=2, nrow=6)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4,	
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=6)

#CARREGA AS MATRIZES NO SIMMR
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))


simmr_out <- simmr_mcmc(simmr_groups)

# Diagnostico de convergencia
summary(simmr_out,type='diagnostics')

#AVALIANDO O AJUSTE DO MODELO
posterior_predictive(simmr_out)

#Comparando as distribuicoes priori e posteriori do modelo Bayesiano
prior_viz(simmr_out, group=4)
# Media e desvios das posterioris
summary(simmr_out,type='statistics', group= 4)

# Quantis das posterioris. 
summary(simmr_out,type='quantiles', group = c(4))

plot(simmr_out,type='density', group = 4, title='Praia das Cabras')

plot(simmr_out,type='boxplot', group = 4, title='Praia das Cabras')

# Probabilidade do ordenamento das contribui??es de cada fonte (na ordem da mais assimilada)
compare_sources(simmr_out,source_names=c("Amarillodesma", "Donax", "Emerita", "Excirolana", "Polichaeta", "Ollivancillaria"), gr= 4)


##### SIMMR para LAGOA DO PEIXE ####
# rodar consumidores do primeiro codigo
#presas
s_names <- c("Amarillodesma", "Donax", "Emerita", "Excirolana", "Polichaeta", "Ollivancillaria", "Tagellus")
s_means <- matrix(c(-14.34499, -14.45716, -13.94067, -13.75714, -12.93064, -12.97285, -14.25975,  
                    10.53560,   10.49058, 10.44184,   11.41262, 10.47439, 13.87949, 10.12610 ), ncol=2, nrow=7) 
s_sds <- matrix(c(0.85304723, 0.64732883, 0.18102714, 0.00557483, 1.07772855, 0.85554136, 0.65922363,
                  0.5478120,  1.1503264,  0.7023311,  0.2082687, 2.1370971, 0.6918654, 1.1388777 ), ncol=2, nrow=7)

#TEF

c_means <- matrix(c(0.2,	0.2,	0.2, 0.2,	0.2, 0.2, 0.2,
                    2.7,	  2.7,  2.7, 2.7,	2.7, 2.7, 2.7), ncol=2, nrow=7)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4, 0.4,	
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=7)
#CARREGA AS MATRIZES NO SIMMR
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))

simmr_out <- simmr_mcmc(simmr_groups)

# Diagnostico de convergencia
summary(simmr_out,type='diagnostics')

#AVALIANDO O AJUSTE DO MODELO
posterior_predictive(simmr_out)

prior_viz(simmr_out, group=5)

summary(simmr_out,type='statistics', group= 5)

summary(simmr_out,type='quantiles', group = 5)

plot(simmr_out,type='density', group = 5, title='Lagoa do Peixe')

plot(simmr_out,type='boxplot', group = 5, title='Lagoa do Peixe')

##### SIMMR MACHOS E FEMEAS PG, grupo1 machos, grupo2 femeas#####

mix <- matrix(c(-14.777, -14.331, -15.802, -13.264, -13.592, -13.742, -13.342,
                -13.859, -13.686, -15.085, -15.292, -15.974, -15.055, -14.918, 
                -15.211, -14.773, -16.383, -13.829, -14.120, -13.726, -13.749, -13.435,
                -13.190, -13.462, -14.009, -13.665, -13.819, -15.748, -15.582, -16.142,
                -14.425, -14.434, 12.352, 12.475, 8.964, 14.192, 13.339, 13.416,
                13.064, 11.800, 12.131, 12.800, 12.324, 12.840, 11.081, 11.220,
                12.289, 14.552, 9.220, 12.750, 12.394, 14.112, 14.308, 13.739, 
                13.600, 13.810, 12.755, 13.185, 12.430, 12.114, 13.212, 11.262,
                14.032, 12.433), ncol=2, nrow=32)

colnames(mix) <-  c("d13C", "d15N")

grp <- as.integer(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))
#carregando dados das presas
s_names <- c("Amarillodesma", "Donax", "Emerita","Excirolana","Olivancillaria", "Perna", "Polichaeta")
s_means <- matrix(c(-14.34499,-14.45716,-13.94067,-13.75714, -12.97285, -14.45223, -12.93064,
                    10.53560, 10.49058, 10.44184, 11.41262, 13.87949, 10.22637, 10.47439), ncol=2, nrow=7) #primeiro carbono depois nitrogenio
s_sds <- matrix(c(0.85304723, 0.64732883,0.18102714, 0.00557483, 0.85554136, 0.73409650, 1.07772855,
                  0.5478120,1.1503264, 0.7023311, 0.2082687, 0.6918654, 1.4103151, 2.1370971 ), ncol=2, nrow=7)

#TEF
c_means <- matrix(c(0.2,	0.2,	0.2,0.2,	0.2,	0.2, 0.2, 
                    2.7,	2.7, 2.7,2.7,	2.7, 2.7, 2.7 ), ncol=2, nrow=7)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4,	0.4,
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=7)

#CARREGA AS MATRIZES NO SIMMR
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))

simmr_out <- simmr_mcmc(simmr_groups)

# Diagnostico de convergencia
summary(simmr_out,type='diagnostics')

#AVALIANDO O AJUSTE DO MODELO
posterior_predictive(simmr_out)

prior_viz(simmr_out, group=1)
prior_viz(simmr_out, group=2)


summary(simmr_out,type='statistics', group= c(1,2))

summary(simmr_out,type='quantiles', group = c(1,2))

plot(simmr_out,type='density', group = 1, title='PG machos')
plot(simmr_out,type='density', group = 2, title='PG femeas')

plot(simmr_out,type='boxplot', group = 1, title='PG machos')
plot(simmr_out,type='boxplot', group = 2, title='PG femeas')


##### SIMMR MACHOS E FEMEAS LP, grupo 1 machos, grupo 2 femeas####

mix <- matrix(c(-13.866535, -13.622131, -14.19175, -14.4607915, -14.258764, -14.0961565, -13.693087,
                -13.8497815, -13.9365055, -14.0409685, -14.138533, -13.894129, 
                11.3317261, 12.388531, 12.7774435, 11.5868527, 11.4292135, 10.9438507, 11.7590113,
                12.3273421, 11.9861362, 12.0856978, 12.0307315, 11.7268612), ncol=2, nrow=12)

colnames(mix) <-  c("d13C", "d15N")

grp <- as.integer(c(1,1,1,1,1,1,
                    2,2,2,2,2,2))
#presas
s_names <- c("Amarillodesma", "Donax", "Emerita", "Excirolana", "Polichaeta", "Ollivancillaria", "Tagellus")
s_means <- matrix(c(-14.34499, -14.45716, -13.94067, -13.75714, -12.93064, -12.97285, -14.25975,  
                    10.53560,   10.49058, 10.44184,   11.41262, 10.47439, 13.87949, 10.12610 ), ncol=2, nrow=7) 
s_sds <- matrix(c(0.85304723, 0.64732883, 0.18102714, 0.00557483, 1.07772855, 0.85554136, 0.65922363,
                  0.5478120,  1.1503264,  0.7023311,  0.2082687, 2.1370971, 0.6918654, 1.1388777 ), ncol=2, nrow=7)

#TEF

c_means <- matrix(c(0.2,	0.2,	0.2, 0.2,	0.2, 0.2, 0.2,
                    2.7,	  2.7,  2.7, 2.7,	2.7, 2.7, 2.7), ncol=2, nrow=7)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4, 0.4,	
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=7)
#CARREGA AS MATRIZES NO SIMMR
simmr_groups <- simmr_load(mixtures=mix,
                           source_names=s_names,
                           source_means=s_means,
                           source_sds=s_sds,
                           correction_means=c_means,
                           correction_sds=c_sds,
                           group=as.factor(paste('consumidor', grp)))

simmr_out <- simmr_mcmc(simmr_groups)

# Diagnostico de convergencia
summary(simmr_out,type='diagnostics')

#AVALIANDO O AJUSTE DO MODELO
posterior_predictive(simmr_out)

prior_viz(simmr_out, group=1)
prior_viz(simmr_out, group=2)


summary(simmr_out,type='statistics', group= 1)
summary(simmr_out,type='statistics', group= 2)

summary(simmr_out,type='quantiles', group = 1)
summary(simmr_out,type='quantiles', group = 2)

plot(simmr_out,type='density', group = 1, title='LP machos')
plot(simmr_out,type='density', group = 2, title='LP femeas')

plot(simmr_out,type='boxplot', group = 1, title='LP machos')
plot(simmr_out,type='boxplot', group = 2, title='LP femeas')

##### SIMMR diferentes ANOS PRAIA GRANDE #####
mix <- matrix(c(-16.3829466, -13.8294, -14.7768832, -14.330768,
                -14.1199888, -15.8021972, -13.7256085, -13.7485672, -13.2635704,
                -13.435468, -13.5920176,-13.1899, -13.4620712, -13.742428,-14.0094832,
                -13.6646648, -13.819168,-13.3423568, -13.8590728, -13.686152, -15.7478545,
                -15.084613, -15.291568, -15.973534, -15.5822905,-16.1420545, -15.055048,
                -14.9180635, -15.210757,-15.5437, -14.4253135, -14.773195, -14.434183,
                -14.7898, -15.1633, -15.2433,-14.6367,-14.9222, -15.0408, -14.6278,
                -14.9963, -14.9341,9.2196328, 12.750148,12.3515168,12.4748048,12.3936402,
                8.9641912,14.1121912,14.3076864,14.1915902,13.7385068,13.3388482,
                13.5998078,13.8104248,13.4159032,12.755285,13.1847382,12.4295992,
                13.063505,11.799803,12.1306258,12.1136995,12.8002597,12.3242308,
                12.8396695,13.2119884,11.2622404,11.0807479,11.2197193,12.2889694,
                13.19394,14.0323345,14.5519216,12.4331263,12.32865,13.19701,12.99654,
                13.49158,13.107,12.99245,10.83741,13.58261,12.972), ncol=2, nrow=42)
colnames(mix) <-  c("d13C", "d15N")

## agrupando os ANOS 1-2018, 2-2019, 3-2020, 4-2021
grp <- as.integer(c(1,1,1,1,1,1,
                    2,2,2,2,2,2,2,
                    2,2,2,2,2,2,2,
                    3,3,3,3,3,3,3,
                    3,3,3,3,3,3,
                    4,4,4,4,4,4,4,
                    4,4))

##### SIMMR diferentes ANOS LAGOA DO PEIXE######

mix <- matrix (c(-13.693087, -13.866535, -14.01436, -13.8497815,
                 -13.9365055, -13.622131, -14.19175, -14.1730255,
                 -14.0409685, -13.8103615, -14.4607915, -14.258764, 
                 -14.138533, -14.0961565, -13.894129, -13.8364, -13.6892,
                 -13.7139, -13.6082, -13.7336, -13.5864, -13.4837,-13.8364,
                 -14.1644, 11.7590113, 11.3317261,11.7818275, 12.3273421, 11.9861362,
                 12.388531, 12.7774435,11.880352,12.0856978,12.0753268,11.5868527,11.4292135,
                 12.0307315,10.9438507, 11.7268612,13.00063,13.3402,13.07632,13.107,
                 12.75925,13.27065,12.97813,12.82471,12.89835
),ncol=2, nrow=24)

grp <- as.integer(c(1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,
                    2,2,2,2,2,2,2,2,2))

#presas
s_names <- c("Amarillodesma", "Donax", "Emerita", "Excirolana", "Polichaeta", "Ollivancillaria", "Tagellus")
s_means <- matrix(c(-14.34499, -14.45716, -13.94067, -13.75714, -12.93064, -12.97285, -14.25975,  
                    10.53560,   10.49058, 10.44184,   11.41262, 10.47439, 13.87949, 10.12610 ), ncol=2, nrow=7) #primeiro carbono depois nitrogenio
s_sds <- matrix(c(0.85304723, 0.64732883, 0.18102714, 0.00557483, 1.07772855, 0.85554136, 0.65922363,
                  0.5478120,  1.1503264,  0.7023311,  0.2082687, 2.1370971, 0.6918654, 1.1388777 ), ncol=2, nrow=7)

#TEF

c_means <- matrix(c(0.2,	0.2,	0.2, 0.2,	0.2, 0.2, 0.2,
                    2.7,	  2.7,  2.7, 2.7,	2.7, 2.7, 2.7), ncol=2, nrow=7)
c_sds <- matrix(c(0.4, 0.4, 0.4, 0.4,	0.4, 0.4, 0.4,	
                  0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), ncol=2, nrow=7)

#### SANKEY DIAGRAM for H. palliatus diet #####
install.packages("networkD3")
library(networkD3)
library(ggplot2)
###### SANKEY bivalves#####
install.packages("networkD3")
library(networkD3)
nodes <- data.frame("name" = 
                      c("Passo de Torres", #Node 0
                        "Praia Grande", # Node 1
                        "Itapeva", # Node 2
                        "Praia das Cabras", # Node 3
                        "Lagoa do Peixe", # Node 4
                        "Emerita", # Node 5
                        "Bivalvia", # Node 6
                        "Excirolana", #Node 7
                        "Polichaeta", #Node 8
                        "Olivancillaria"#Node 9
                      )) 
#
links <- as.data.frame(matrix(c(
  0, 6, 36.7,  
  0, 5, 17.4,
  0, 7, 13,
  0, 9, 8.5,
  0, 8, 24.6,
  1, 6, 88.6,  
  1, 5, 4.2,
  1, 7, 2.9,
  1, 9, 1.7,
  1, 8, 2.6,
  2, 6, 58.7,  
  2, 5, 10.4,
  2, 7, 14.3,
  2, 9, 10.9,
  2, 8, 5.7,
  3, 6, 60.6,
  3, 5, 17.2,
  3, 7, 9.6,
  3, 9, 4.9,
  3, 8, 7.7,
  4, 6, 53.7,
  4, 5, 31.4,
  4, 7, 6.5,
  4, 9, 2.5,
  4, 8, 5.9 
), byrow = TRUE, ncol = 3))

names(links) <- c("source", "target", "value")
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30)

#### SANKEY por especie #####
nodes <- data.frame("name" = 
                      c("Passo de Torres", #Node 0
                        "Praia Grande", # Node 1
                        "Itapeva", # Node 2
                        "Praia das Cabras", # Node 3
                        "Lagoa do Peixe", # Node 4
                        "Emerita", # Node 5
                        "Donax", # Node 6
                        "Amarilladesma", # Node 7
                        "Excirolana", #Node 8
                        "Polichaeta", #Node 9
                        "Olivancillaria",#Node 10
                        "Perna", #Node 11
                        "Tagellus" #Node 12
                      )) 

links = as.data.frame(matrix(c(
  0, 7, 12.7, # Each row represents a link. The first number
  0, 6, 11.8, # represents the node being conntected from. 
  0, 5, 17.4,
  0, 8, 13,
  0, 10, 8.4,
  0, 11, 12.2,
  0, 9, 24.6,
  1, 7, 6.7, 
  1, 6, 13.2,  
  1, 5, 4.2,
  1, 8, 2.9,
  1, 10, 1.7,
  1, 11, 68.7,
  1, 9, 2.6,
  2, 7, 21.9,
  2, 6, 21.5, 
  2, 5, 10.4,
  2, 8, 14.3,
  2, 10, 10.9,
  2, 11, 15.3,
  2, 9, 5.7,
  3, 7, 24.6,
  3, 6, 36,
  3, 5, 17.2,
  3, 8, 9.6,
  3, 9, 7.7,
  3, 10, 4.9,
  4, 7, 11.9,
  4, 6, 17.2,
  4, 5, 31.4,
  4, 9, 5.9,
  4, 12, 24.6,
  4, 8, 6.5,
  4, 10, 2.5
), # the second number represents the node connected to.
byrow = TRUE, ncol = 3))

names(links) = c("source", "target", "value")
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 20)
#### PREY: ABUNDANCE, NMDS, PERMANOVA ####
dados <- read.csv(file = "macrobentos_abundancia_primok.csv", header=T, sep=";") #"macrobentos_abundancia_primok.csv"
head(dados) 
dados[dados > 0] <- 1 
dados
########## nMDS ########### 

library (vegan)

initial_nMDS1 <- metaMDS(dados[,3:8], distance="bray", k=2, trymax=1000)
nMDS1 <- metaMDS(dados[,3:8], previous.best = initial_nMDS1, k=2, trymax=1000)

stressplot(nMDS1)

#nMDS x local
op <- par(mar=c(4,4,1,1))
plot(nMDS1$points[,1], nMDS1$points[,2], font.lab=2, pch=16, 
     bg=c("black", "orange")
     [dados$Local], xlab='NMDS1', ylab='NMDS2', cex=0.9, 
     xlim=c(-2, 2), ylim=c(-2, 2))

ordihull(nMDS1,groups=dados$Local,draw="polygon",col=c("green", "yellow", "magenta", "orange"),
         label=F, air=2)


ordiellipse(nMDS1,groups=dados$Local,draw="polygon",col=c("magenta", "orange", "green", "yellow"),
            label=T, conf=0.7, kind="sd", alpha = 30)
ordilabel(x=nMDS1$species, dis="sites", cex=0.8, font=3, fill="white", col="blue")

#### ABUNDANCIA ####

library(vegan)
dados <- read.csv(file = "macrobentos_abundancia_primok.csv", header=T, sep=";")

attach(dados)

pg1 <- dados [which(dados$Local=="PG19"),]
pg2 <- dados [which(dados$Local=="PG20"),]
peva1 <- dados[which(dados$Local=="ITA19"),]
peva2 <- dados[which(dados$Local=="ITA20"),]
cb <- dados[which(dados$Local=="CB19"),]
pnlp1 <- dados[which(dados$Local=="LP19"),]
pnlp2 <- dados[which(dados$Local=="LP20"),]

#somar a abundancia de todos os pontos por especie

abund.pg1 <- colSums(pg1 [,3:8])
abund.pg2 <- colSums(pg2 [,3:8])
abund.peva1 <- colSums(peva1 [,3:8])
abund.peva2 <- colSums(peva2 [,3:8])
abund.cb <- colSums(cb [,3:8])
abund.pnlp1 <- colSums(pnlp1 [,3:8])
abund.pnlp2 <- colSums(pnlp2 [,3:8])

abund.primavera <- cbind(abund.pg1, abund.pg2, abund.peva1, 
                         abund.peva2, abund.cb, abund.pnlp1, abund.pnlp2)

abund.primavera <- as.matrix(abund.primavera)
abund.primavera <- t(abund.primavera)
abundancia <- abund.primavera
abundancia <- t(abundancia)
abundancia <- as.data.frame(abundancia)
# calculando a abundancia relativa por local #

relativa.pg1 <- (abundancia$abund.pg1*100/sum(abundancia$abund.pg1))
relativa.pg2 <- (abundancia$abund.pg2*100/sum(abundancia$abund.pg2))
relativa.peva1 <- (abundancia$abund.peva1*100/sum(abundancia$abund.peva1))
relativa.peva2 <- (abundancia$abund.peva2*100/sum(abundancia$abund.peva2))
relativa.cb <- (abundancia$abund.cb*100/sum(abundancia$abund.cb))
relativa.pnlp1 <- (abundancia$abund.pnlp1*100/sum(abundancia$abund.pnlp1))
relativa.pnlp2 <- (abundancia$abund.pnlp2*100/sum(abundancia$abund.pnlp2))

abundancia.relativa <- cbind(relativa.pg1, relativa.pg2, relativa.peva1, relativa.peva2, 
                             relativa.cb, relativa.pnlp1, relativa.pnlp2)

tabela.abundancia <- cbind(abundancia, abundancia.relativa)
tabela.abundancia.t <- as.data.frame(t(tabela.abundancia))


# carregar arquivo csv de resultados # 
abund.relativa <- read.csv(file="abund.rel.primavera.csv", head=T, sep=";")

attach(abund.relativa)

#### barplot relative abundance####
abund.relativa.matrix <- as.matrix(abund.relativa)
par(mar=c(9,4,2,2))
barplot(abund.relativa.matrix[,-1], col= c("paleturquoise1", "rosybrown1",
                                           "palegreen", "skyblue",
                                           "navajowhite", "orange"))

legend(x=2, y=-14, xpd=TRUE, ncol=2, 
       legend=c("E_brasiliensis", "D_hanleyanus", "A_mactroides","Amphipoda",
                "E_armata", "Polychaeta"),
       fill=c("paleturquoise1", "rosybrown1",
              "palegreen", "skyblue",
              "navajowhite", "orange"), bty="n")
#### FREQUENCIA DE OCORRENCIA ##### 
# transformar os dados em presenca e ausencia#
# calculando a frequencia de ocorrencia de cada especie 
library(plyr)
dados <- read.csv(file = "macrobentos_abundancia_primok.csv", header=T, sep=";") #"macrobentos_abundancia_primok.csv"
freq.oc.primavera <- ldply(dados[,3:8], function(c) sum(c=="1"))
fo.relativa <- freq.oc.primavera[,2]*(100/70)
# ocorrencia na praia grande
fo.pg1 <- ldply(dados[1:10,3:8], function(c) sum(c=="1"))
fo.rel.pg1 <- (fo.pg1[,2]*100)/sum(fo.pg1[,2])
fo.pg2 <- ldply(dados[11:20,3:8], function(c) sum(c=="1"))
fo.rel.pg2 <- (fo.pg2[,2]*100)/sum(fo.pg2[,2])
# ocorrencia itapeva
fo.peva1 <- ldply(dados[21:30,3:8], function(c) sum(c=="1"))
fo.rel.peva1 <- (fo.peva1[,2]*100)/sum(fo.peva1[,2])
fo.peva2 <- ldply(dados[31:40,3:8], function(c) sum(c=="1"))
fo.rel.peva2 <- (fo.peva2[,2]*100)/sum(fo.peva2[,2])
# ocorrencia praia das cabras
fo.cb <- ldply(dados[41:50,3:8], function(c) sum(c=="1"))
fo.cb.rel <- fo.cb[,2]*(100/10)#fo relativa
fo.rel.cb <- (fo.cb[,2]*100)/sum(fo.cb[,2])
#ocorrencia lagoa do peixe
fo.pnlp1 <- ldply(dados[51:60,3:8], function(c) sum(c=="1"))
fo.rel.pnlp1 <- (fo.pnlp1[,2]*100)/sum(fo.pnlp1[,2])
fo.pnlp2 <- ldply(dados[61:70,3:8], function(c) sum(c=="1"))
fo.rel.pnlp2 <- (fo.pnlp2[,2]*100)/sum(fo.pnlp2[,2])

#### grafico com os dados de frequencia de ocorrencia ####
fo.separados <- cbind(fo.pg1[,1], fo.rel.pg1, fo.rel.peva1, fo.rel.cb, fo.rel.pnlp1, 
                      fo.rel.pg2, fo.rel.peva2, fo.rel.pnlp2)
colnames(fo.separados) <- c("species", "pg1", "peva1", "cb", "pnlp1", "pg2", "peva2", "pnlp2")
fo.separados <- as.matrix(fo.separados)

freq.ocorrencia <- read.csv(file="freq.oc.primavera.csv", head=T, sep=";")
freq.ocorrencia <- as.matrix(freq.ocorrencia)


par(mar=c(9,4,2,2))
barplot(freq.ocorrencia[,-1], col= c("paleturquoise1", "rosybrown1",
                                     "palegreen", "skyblue",
                                     "navajowhite", "orange"))

legend(x=2, y=-12, xpd=TRUE, ncol=2, 
       legend=c("E_brasiliensis", "D_hanleyanus", "A_mactroides","Amphipoda",
                "E_armata", "Polychaeta"),
       fill=c("paleturquoise1", "rosybrown1",
              "palegreen", "skyblue",
              "navajowhite", "orange"), bty="n")

#### PERMANOVA ####
dados <- read.csv(file = "macrobentos_abundancia_primok.csv", header=T, sep=";") #"macrobentos_abundancia_primok.csv"

permanova <- adonis(dados[,3:8] ~ dados$Local, permutations=999, 
                    distance='bray')
# post hoc para permanova
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", force=TRUE) 
install.packages("pairwiseAdonis")
library(pairwiseAdonis)

pair.mod<-pairwise.adonis(dados[,3:8],factors=dados$Local)
pair.mod

# permanova 2019 #
dados <- read.csv(file = "mbentosabundancia_2019.csv", header=T, sep=";") 

permanova2019 <- adonis(dados[,3:8] ~ dados$Local, permutations=999, 
                        distance='bray')

pair.mod<-pairwise.adonis(dados[,3:8],factors=dados$Local)

# permanova 2020 # 
dados <- read.csv(file = "mbentosabundancia_2020.csv", header=T, sep=";") 

permanova2020 <- adonis(dados[,3:8] ~ dados$Local, permutations=999, 
                        distance='bray')

pair.mod<-pairwise.adonis(dados[,3:8],factors=dados$Local)

#### nmds 2019 #####
dados <- read.csv(file = "mbentosabundancia_2019.csv", header=T, sep=";") #"macrobentos_abundancia_primok.csv"
head(dados) 

library (vegan)

initial_nMDS1 <- metaMDS(dados[,3:8], distance="bray", k=2, trymax=1000) #selecionar colunas das especies
nMDS1 <- metaMDS(dados[,3:8], previous.best = initial_nMDS1, k=2, trymax=1000)

stressplot(nMDS1)

op <- par(mar=c(4,4,1,1))
plot(nMDS1$points[,1], nMDS1$points[,2], font.lab=2, pch=16, 
     bg=c("black", "orange")
     [dados$Local], xlab='NMDS1', ylab='NMDS2', cex=0.9, 
     xlim=c(-2, 2), ylim=c(-2, 2))

ordiellipse(nMDS1,groups=dados$Local,draw="polygon",col=c("magenta", "orange", "green", "yellow"),
            label=T, conf=0.7, kind="sd", alpha = 30)

###### nmds 2020
dados <- read.csv(file = "mbentosabundancia_2020.csv", header=T, sep=";") #"macrobentos_abundancia_primok.csv"
head(dados) 

library (vegan)

initial_nMDS1 <- metaMDS(dados[,3:8], distance="bray", k=2, trymax=1000) #selecionar colunas das especies
nMDS1 <- metaMDS(dados[,3:8], previous.best = initial_nMDS1, k=2, trymax=1000)

stressplot(nMDS1)

op <- par(mar=c(4,4,1,1))
plot(nMDS1$points[,1], nMDS1$points[,2], font.lab=2, pch=16, 
     bg=c("black", "orange")
     [dados$Local], xlab='NMDS1', ylab='NMDS2', cex=0.9, 
     xlim=c(-1, 1), ylim=c(-1, 1))

ordiellipse(nMDS1,groups=dados$Local,draw="polygon",col=c("magenta", "orange", "green", "yellow"),
            label=T, conf=0.7, kind="sd", alpha = 30)
