
library(ggplot2)
library(egg)
library(minpack.lm)
library(brms)
library(propagate)
library(broom)
library(broom.mixed)
library(dplyr)
library(plyr)
library(ggdistribute)
library(viridis)
library(patchwork)
library(cowplot)
library(ggnewscale)
library(stringr)
library(ggstance)
library(bobfunctions2)


# Dry-wet weight conversion coefficients
conv_dat <- read.csv("dry-wet_conversion_coefficients.csv", sep = ",")




# Load buoyant weight data
dat <- read.csv("Buoyant_weight_Kelso_Reef.csv", sep = ",")
colnames(dat)[1] <- "Plug_number"

n <- grep("Date", colnames(dat))

dat[ ,n] <- lapply(dat[ , n], as.Date, format = "%d/%m/%Y")



# Get data in the correct format (i.e., one column with all buoyant weights and
# one column with the date)
meta <- data.frame("Plug" = dat$Plug_number,
                   "Species" = dat$Species,
                   "position" = dat$position,
                   "tank" = dat$tank,
                   "Dry_weight" = dat$Dry.weight,
                   "Colony" = dat$Colony,
                   "Temperature" = dat$Temperature,
                   "Reef" = dat$Reef,
                   "Notes_fragging" = dat$Notes_fragging)

meta$Buoyant_weight <- NA
meta$Date <- NA
meta$Notes <- NA
meta$Time <- NA


dat2 <- meta

for (i in 1 : nrow(meta)) {
  sub_d <- dat[i, ]
  dat4 <- meta[1, ]
  
  for (j in n) {
    
    dat3 <- meta[i ,]
    dat3$Buoyant_weight <- sub_d[ , j-1]
    dat3$Date <- as.character(sub_d[, j])
    dat3$Notes <- sub_d[ ,j+1]
    
    dat4 <- rbind(dat4, dat3)
    
  }
  
  dat4 <- dat4[-1, ]
  dat4 <- dat4[!is.na(dat4$Buoyant_weight) == TRUE, ]
  dat4$Time <- c(1:nrow(dat4))
  
  dat2 <- rbind(dat2, dat4)
}

dat2 <- dat2[ - c(1:nrow(meta)), ]

# Get survival
dat2$Alive <- 1
dat2$Alive[grep("Sampled", dat2$Notes, ignore.case = TRUE)] <- 0
dat2$Alive[grep("Not Sampled", dat2$Notes, ignore.case = TRUE)] <- 1
dat2[dat2$position == 0, ]$Alive <- NA

# Label plugs (always in position 0) as plugs
dat2[dat2$position == 0, ]$Species <- "plug"
dat2[dat2$position == 0, ]$Reef <- "plug"

dat2 <- dat2[complete.cases(dat2$Buoyant_weight), ]
dat2 <- dat2[!dat2$Buoyant_weight == "", ]


## Days relative to first day sampled
dat2$Days <- NA
for (i in 1:nrow(dat2)) {
  
  sub_d <- dat2[dat2$Plug == dat2$Plug[i], ]
  dat2$Days[i] <- as.numeric(difftime(dat2$Date[i], min(sub_d$Date), 
                                      units = "days"))
}

dat2$Buoyant_weight <- as.numeric(dat2$Buoyant_weight)
dat2$Buoyant_weight <- round(dat2$Buoyant_weight, digits = 2)



# Create a column to identify data with problems
dat2$Problem <- 0

# After data exploration, there are obvious typos in :
dat2[dat2$Plug == 482 & dat2$Time == 1, ]$Problem <- "typo"
dat2[dat2$Plug == 241 & dat2$Time == 1, ]$Problem <- "typo" 
dat2[dat2$Plug == 278 & dat2$Time == 1, ]$Problem <- "typo" 
dat2[dat2$Plug == 5   & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 53  & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 365 & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 219 & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 393 & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 299 & dat2$Time == 5, ]$Problem <- "typo"
dat2[dat2$Plug == 199 & dat2$Time == 1, ]$Problem <- "typo"
dat2[dat2$Plug == 63  & dat2$Time == 1, ]$Problem <- "typo"
dat2[dat2$Plug == 409 & dat2$Time == 5, ]$Problem <- "dead"
dat2[dat2$Plug == 67  & dat2$Time == 3, ]$Problem <- "dead"


# Eliminate observations that have problems with excess algae or breakage of 
# branches
dat2[dat2$Plug == 466 & dat2$Time == 3, ]$Problem <- "Algae"
dat2[dat2$Plug == 416 & dat2$Time == 3, ]$Problem <- "Algae"
dat2[dat2$Plug == 273 & dat2$Time == 2, ]$Problem <- "Algae"
dat2[dat2$Plug == 528 & dat2$Time == 3, ]$Problem <- "Algae"
dat2[dat2$Plug == 511 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 263 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 315 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 129 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 41  & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 608 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 199 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 568 & dat2$Time == 5, ]$Problem <- "Algae"
dat2[dat2$Plug == 193 & dat2$Time  > 5, ]$Problem <- "Breakage"
dat2[dat2$Plug == 458 & dat2$Time  > 3, ]$Problem <- "Breakage"
dat2[dat2$Plug == 520 & dat2$Time  > 2, ]$Problem <- "Breakage"


# Exclude problematic data points
Dat <- dat2[dat2$Problem == 0, ]



# Function to convert plug dry weight to wet weight
convert_to_wet <- function(dry_weight, temp, reg_coef) {
  
  reg_ <- reg_coef[reg_coef$Temperature == temp, ]
  
  wet_weight <- reg_$intercept + reg_$slope * dry_weight
  
  return(wet_weight)
  
} 

Dat$Wet_weight_plug <- NA

# Convert plug dry weight to wet weight
for (i in 1:nrow(Dat)) {
  Dat$Wet_weight_plug[i] <- convert_to_wet(Dat$Dry_weight[i], 
                                           Dat$Temperature[i], conv_dat)
}



# Dry weight of plug 564 must be wrong
Dat[Dat$Wet_weight_plug < 6, ]$Wet_weight_plug <- NA

# Function to get fragment dry weight
convert_to_dry <- function(Species, wet_weight_plug, buoyant_weight, temp) {
  
  # seawater density (salinity is 35psu)
  skeletal_density <- ifelse(
    Species == "Acropora hyacinthus", 1.29, # Mohammed & Dar (2017)
    ifelse(
      Species == "Acropora tenuis", 2.55,     # Nielsen et al (2020)
      ifelse(
        Species == "Pocillopora verrucosa", 2.78, # Spencer Davis (1989)
        1.62 # Stylophora pistillata: Mohammed & Dar (2017)
      )
      
    ))
  seawater_dens <- data.frame("Temp" = c(19,21,23,25,26,27,28,29,30,31),
                              "Dens" = c(1.029,1.029,1.028,1.028,1.027,
                                         1.027,1.027,1.026,1.026,1.026))
  density <- seawater_dens[which(seawater_dens$Temp == temp), ]$Dens
  
  dry_weight = (buoyant_weight - wet_weight_plug) /
    (1 - density/skeletal_density)
  
  if (Species == "plug") {
    dry_weight = (buoyant_weight) /
      (1 - density/2.95)
    
  }
  
  return(dry_weight)
}


# Get growth
growth <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(growth) <- c("Plug", "Species", "Reef", "Colony", "Temperature", "Tank", 
                      "Position", "Dry_weight_plug", "Wet_weight_plug",
                      "Initial_weight", "Initial_day", "Initial_notes", 
                      "Final_weight", "Final_day", "Final_notes", "Max_span")




for (p in unique(Dat$Plug)) {
  
  sub_d <- Dat[Dat$Plug == p, ]
  
  if (nrow(sub_d) > 1) {
    
    for (t in 1:(nrow(sub_d) - 1)) {
      for (t1 in (t+1): nrow(sub_d)) {
        
        growth1 <- data.frame(
          "Plug" = sub_d$Plug[t],
          "Species" = sub_d$Species[t],
          "Reef" = sub_d$Reef[t],
          "Colony" = sub_d$Colony[t],
          "Temperature" = sub_d$Temperature[t],
          "Tank" = sub_d$tank[t],
          "Position" = sub_d$position[t],
          "Dry_weight_plug" = sub_d$Dry_weight[t],
          "Wet_weight_plug" = sub_d$Wet_weight_plug[t],
          "Initial_weight" = sub_d$Buoyant_weight[t],
          "Initial_day" = sub_d$Days[t],
          "Initial_notes" = sub_d$Notes[t], 
          "Final_weight" = sub_d$Buoyant_weight[t1],
          "Final_day" = sub_d$Days[t1],
          "Final_notes" = sub_d$Notes[t1],
          "Max_span" =  ifelse(t == 1 & t1 == nrow(sub_d), "Y", "N")
        )
        
        growth <- rbind(growth, growth1)
      }
    }
  }
  
  
  
}


# Calculate change in dry weight per day in g
for (i in 1:nrow(growth)) {
  
  if (growth$Species[i] == "plug") {
    growth$Initial_weight[i] == growth$Dry_weight_plug[i]
    growth$Final_weight[i] == growth$Final_weight[i]
  }
  
  else{
    growth$Initial_weight[i] <- convert_to_dry(growth$Species[i],
                                               growth$Wet_weight_plug[i],
                                               growth$Initial_weight[i],
                                               growth$Temperature[i])
    growth$Final_weight[i] <- convert_to_dry(growth$Species[i],
                                             growth$Wet_weight_plug[i],
                                             growth$Final_weight[i],
                                             growth$Temperature[i])
  }
  
}





growth$Time <- growth$Final_day - growth$Initial_day

Growth <- growth[growth$Time > 12 & growth$Time < 17, ]
Growth$mid_point <- (Growth$Final_day + Growth$Initial_day)/2
Growth$mid_week <- round(Growth$mid_point/7, 0)


Growth$Growth <- (Growth$Final_weight - Growth$Initial_weight) / Growth$Time


# Get growth relative to initial weight
Growth$Growth <- Growth$Growth / Growth$Initial_weight

# Get relative growth in mg instead of g
Growth$Growth <- Growth$Growth * 1000

Growth$mid_point <- (Growth$Final_day + Growth$Initial_day)/2
Growth$mid_week <- round(Growth$mid_point/7, 0)


Growth <- Growth[!Growth$Species == "plug", ]

############## Time effect #########################
colour_values <-  c("#999999", "#E69F00", "#56B4E9", "#009E73")
names(colour_values) <- c("Acropora hyacinthus",
                          "Acropora tenuis",
                          "Pocillopora verrucosa",
                          "Stylophora pistillata")


Growth$colony <- paste(Growth$Species, Growth$Colony)
Growth$Species <- factor(Growth$Species)

mod0_2 <- readRDS("TPC_2022_0_2.rds")

mod2_4 <- brm(bf(Growth ~ Gmax * exp( -exp(ro * (Temperature - xshift)) - exp(sigm) * (Temperature - xshift)^2),
              Gmax ~  Species - 1 + (1 | Species : colony),
              ro   ~  Species - 1 + (1 | Species : colony),
              xshift ~  Species - 1 + (1 | Species : colony),
              sigm ~  Species - 1 + (1 | Species : colony), nl = TRUE), 
           data = Growth[Growth$Reef == "Kelso" & Growth$mid_week == 3, ],
           family = gaussian(),
           prior = c(
             prior(normal(0, 10), nlpar = Gmax),
             prior(normal(0, 5), nlpar = ro),
             prior(normal(25, 10), nlpar = xshift),
             prior(normal(0, 5), nlpar = sigm)),
           control = list(max_treedepth = 18, adapt_delta = 0.99),
           chains = 3, cores = 3, iter = 10e4, 
           file = "TPC_2022_2_4")


mod4_6 <- brm(bf(Growth ~ Gmax * exp( -exp(ro * (Temperature - xshift)) - exp(sigm) * (Temperature - xshift)^2),
                 Gmax ~  Species - 1 + (1 | Species : colony),
                 ro   ~  Species - 1 + (1 | Species : colony),
                 xshift ~  Species - 1 + (1 | Species : colony),
                 sigm ~  Species - 1 + (1 | Species : colony), nl = TRUE), 
              data = Growth[Growth$Reef == "Kelso" & Growth$mid_week == 5, ],
              family = gaussian(),
              prior = c(
                prior(normal(0, 10), nlpar = Gmax),
                prior(normal(0, 5), nlpar = ro),
                prior(normal(25, 10), nlpar = xshift),
                prior(normal(0, 5), nlpar = sigm)),
              control = list(max_treedepth = 18, adapt_delta = 0.99),
              chains = 3, cores = 3, iter = 10e4, 
              file = "TPC_2022_4_6")


df <- expand.grid("Species" = unique(Growth$Species), 
                  "Temperature" = seq(19, 31, by = 0.1))

fitted0_2 <- cbind(df, as.data.frame(fitted(mod0_2, df, re_formula = NA)))
fitted0_2$mid_week <- 1
fitted2_4 <- cbind(df, as.data.frame(fitted(mod2_4, df, re_formula = NA)))
fitted2_4$mid_week <- 3
fitted4_6 <- cbind(df, as.data.frame(fitted(mod4_6, df, re_formula = NA)))
fitted4_6$mid_week <- 5

fitted <- rbind(fitted0_2, fitted2_4, fitted4_6)

p <- ggplot() +
  geom_jitter(data = Growth[Growth$Reef == "Kelso" & Growth$mid_week %in% c(1, 3, 5), ],
             aes( x = Temperature, y = Growth, group = mid_week, fill = mid_week),
             pch = 21, col = "black", height = 0, width = 0.2) +
  geom_ribbon(data = fitted, aes( x = Temperature, ymin = Q2.5, ymax = Q97.5,
                                  fill = mid_week, group = mid_week),
              alpha = 0.3) +
  geom_line(data = fitted, aes( x = Temperature, y = Estimate,
                                  col = mid_week, group = mid_week),
              size = 1.1, show.legend = FALSE) +
  facet_wrap( ~ Species) + 
  theme_bw() +
  scale_colour_viridis(breaks = c(1, 3, 5), name = "Weeks",
                       labels = c("0-2", "2-4", "4-6")) + 
  scale_fill_viridis(breaks = c(1, 3, 5), name = "Weeks",
                       labels = c(" 0-2", "2-4", "4-6")) + 
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  geom_hline(yintercept = 0, col = "black") +
  scale_x_continuous( breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32),
                      labels = c(19, 21, 23, 25, "", 27, "", 29, "", 31, "") ) +
  xlab( "Temperature (°C)") +
  theme(strip.text = element_text(face = "italic")) +
  
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text.align = 0) 
p
                   