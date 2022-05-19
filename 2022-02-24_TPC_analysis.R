setwd( "~/Temperature-Growth/Data")

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
conv_dat <- read.csv("dry-wet_conversion_coeffs.csv", sep = ",")

# Load temperature data at the site from AIMS data loggers
temp <- read.csv("Kelso_water_temp_pday.csv", sep = ",")
temp <- temp[!temp$date < "1996-01-01", ]
temp <- temp[!temp$date > "2019-12-31", ]
temp$year <- c()
temp$month <- c()
temp$day <- c()

for( i in 1:nrow(temp)) {
  temp$year[i] <- unlist(strsplit(temp$date[i], "-"))[1]
  temp$month[i] <- unlist(strsplit(temp$date[i], "-"))[2]
  temp$day[i] <- unlist(strsplit(temp$date[i], "-"))[3]
}

temp <- temp[complete.cases(temp$depth8m), ]

# Calculate number of temperature daily records per year
temp1 <- temp %>% 
  dplyr::group_by(year) %>%
  dplyr::summarise("n" = n())
temp1 <- as.data.frame(temp1)

temp <- merge(temp, temp1, by = "year" , all.x = TRUE)

# Exclude years with less than 360 days of data
sub_temp <- temp[temp$n >= 360, ]
sub_temp$days <- as.Date(sub_temp$date, "%Y-%m-%d") - as.Date(min(sub_temp$date),"%Y-%m-%d")
sub_temp$month_day <- paste(sub_temp$month, sub_temp$day, sep = "-")


# Get median temperature for each day of the year
sum_temp <- sub_temp %>% dplyr::group_by(month_day) %>% 
  dplyr::summarise("median" = median(depth8m, na.rm =TRUE),
                   "upr" = quantile(depth8m, 0.75, na.rm = TRUE),
                   "lwr" = quantile(depth8m, 0.25, na.rm = TRUE))
sum_temp <- as.data.frame(sum_temp)



day_num <- data.frame("month_day" = unique(sum_temp$month_day))
day_num$number_day <- c(1:366)
sum_temp <- merge(sum_temp, day_num, by = "month_day")
sub_temp <- merge(sub_temp,  day_num, by = "month_day")

# Plot temperature data across the year. Red shows summer days, 
# blue shows winter days
FigS7 <- 
  ggplot() +
  geom_rect(aes(xmin = 336, xmax= 366, ymin = -Inf, ymax = Inf), fill = "red", 
            alpha = 0.2) +
  geom_rect(aes(xmin = 1, xmax= 60, ymin = -Inf, ymax = Inf), fill = "red", 
            alpha = 0.2) +
  geom_rect(aes(xmin = 153, xmax= 244, ymin = -Inf, ymax = Inf), fill = "blue", 
            alpha = 0.2) +
  geom_point(data = sub_temp, aes(x = number_day, y = depth8m), pch = 21, 
             fill = "grey", col = "black") +
  geom_line(data = sum_temp, aes(x = number_day, y = median), size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306,
                                336),
                     labels = c("1-Jan", "1-Feb", "1-Mar", "1-Apr",
                                "1-May", "1-Jun", "1-Jul", "1-Aug",
                                "1-Sep", "1-Oct", "1-Nov", "1-Dec")) +
  ylab("Water temperature at Kelso Reef (°C) (8m-depth)") +
  xlab("") +
  theme(text = element_text(size = 15))
FigS7


# Load buoyant weight data
dat <- read.csv("Buoyant_weight_final_data.csv", sep = ",")
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
dat2[dat2$Plug == 63 & dat2$Time == 1, ]$Problem <- "typo"
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



# Select only growth in the first four weeks
Growth <- growth[growth$Final_day < 30, ]

d <- Growth[1, ]
d[1, ] <- NA

# Find growth estimates using the latest data point (~4 weeks)
for (plug in unique(Growth$Plug)) {
  
  sub_d <- Growth[Growth$Plug == plug, ]
  
  if (nrow(sub_d) == 1 ) {
    d <- rbind(d, sub_d)
  } else {
    n <- which(sub_d$Time == max(sub_d$Time))
    d <- rbind(d, sub_d[n, ])
  }
}

d <- d[-1, ]



# Assign time between measurements as 27.5 days for those fragments that 
# died before four weeks
d$Growth <- ifelse(d$Time < 26, (d$Final_weight - d$Initial_weight) / 27.5,
                   (d$Final_weight - d$Initial_weight) / d$Time)


# Get growth relative to initial weight
d$Growth <- d$Growth / d$Initial_weight

# Get relative growth in mg instead of g
d$Growth <- d$Growth * 1000

plugs <- d[d$Species == "plug", ]

ggplot(plugs) +
  geom_boxplot(aes(y = Growth), fill = "grey") +
  geom_jitter(aes(x = 0, y = Growth), width = 0.2, height = 0,
              pch = 21, fill = "red", col = "black", size = 5, alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = c(), limits = c(-1, 1)) +
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  xlab( "") +
  
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text.align = 0)


d_Kelso <- d[d$Reef == "Kelso", ]
#write.table(d_Kelso, "Coral_fragment_growth_Kelso_Reef_2021.csv", sep =",", row.names = FALSE)

d <- d[!d$Species == "plug" , ]



FigS6 <- ggplot(d) +
  geom_histogram(aes(x = Time), col = "black", fill = "grey") +
  theme_bw() +
  xlab("Days between initial and final measurement") +
  ylab("Number of fragments") +
  scale_x_continuous(breaks = seq(0, 30, by = 5))+
  geom_vline(xintercept = 27.5, col = "red", linetype = "dotted", size = 1.5) +
  theme(text = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0, 300, by = 50))
FigS6




# Get colony id
d$colony <- paste(d$Species, d$Colony)
d$Species <- factor(d$Species)





mod <- brm(bf(Growth ~ Gmax * exp( -exp(ro * (Temperature - xshift)) - exp(sigm) * (Temperature - xshift)^2),
              Gmax ~  Species - 1 + (1 | Species : colony),
              ro   ~  Species - 1 + (1 | Species : colony),
              xshift ~  Species - 1 + (1 | Species : colony),
              sigm ~  Species - 1 + (1 | Species : colony), nl = TRUE), 
           data = d[d$Reef == "Kelso", ],
           family = gaussian(),
           prior = c(
             prior(normal(0, 10), nlpar = Gmax),
             prior(normal(0, 5), nlpar = ro),
             prior(normal(25, 10), nlpar = xshift),
             prior(normal(0, 5), nlpar = sigm)),
           control = list(max_treedepth = 17, adapt_delta = 0.99),
           chains = 3, cores = 3, iter = 5e4, thin = 5, 
           file = "TPC_2022_0_2")

# Create new data frame to make predictions
d_Kelso <- expand.grid("Temperature" = seq(19, 31, by = 0.1),
                       "Species" = levels(d$Species))

d_Kelso <- cbind(d_Kelso, as.data.frame(fitted(mod, d_Kelso, re_formula = NA)))


colours <- c("#999999", "#E69F00", "#56B4E9", "#009E73" )
names(colours) <- c("Acropora hyacinthus", "Acropora tenuis",
                    "Pocillopora verrucosa", "Stylophora pistillata")

# Get proportion of days that fell within each temperature interval
sub_temp$interval <- ifelse(sub_temp$depth8m < 22, "[21-22)",
                            ifelse(sub_temp$depth8m < 23, "[22-23)",
                                   ifelse(sub_temp$depth8m < 24, "[23-24)",
                                          ifelse(sub_temp$depth8m < 25, "[24-25)",
                                                 ifelse(sub_temp$depth8m < 26, "[25-26)",
                                                        ifelse(sub_temp$depth8m < 27, "[26-27)",
                                                               ifelse(sub_temp$depth8m < 28, "[27-28)",
                                                                      ifelse(sub_temp$depth8m < 29, "[28-29)",
                                                                             ifelse(sub_temp$depth8m < 30, "[29-30)",
                                                                                    ifelse(sub_temp$depth8m < 31, "[30-31)",
                                                                                           "[31-32)"))))))))))
sub_temp$interval_midpoint <- ifelse(sub_temp$depth8m < 22, 21.5,
                                     ifelse(sub_temp$depth8m < 23, 22.5,
                                            ifelse(sub_temp$depth8m < 24, 23.5,
                                                   ifelse(sub_temp$depth8m < 25, 24.5,
                                                          ifelse(sub_temp$depth8m < 26, 25.5,
                                                                 ifelse(sub_temp$depth8m < 27, 26.5,
                                                                        ifelse(sub_temp$depth8m < 28, 27.5,
                                                                               ifelse(sub_temp$depth8m < 29, 28.5,
                                                                                      ifelse(sub_temp$depth8m < 30, 29.5,
                                                                                             ifelse(sub_temp$depth8m < 31, 30.5,
                                                                                                    31.5))))))))))
intervals <- as.data.frame(table(sub_temp$interval_midpoint))
intervals$temp <- as.numeric(as.character(intervals$Var1))
intervals$y_val <- -0.13
intervals$proportion <- intervals$Freq / sum(intervals$Freq)
intervals$rel_freq <- intervals$Freq / max(intervals$Freq)
intervals$num_days <- round(intervals$proportion * 365, 0)


Fig1 <-
  ggplot() +
  geom_ribbon(data = d_Kelso, aes(x = Temperature, ymin = Q2.5, ymax = Q97.5, 
                                  fill = Species), alpha = 0.2) +
  geom_line(data = d_Kelso, aes(x = Temperature, y = Q2.5, col = Species), 
            linetype = "dotted") +
  geom_line(data = d_Kelso, aes(x = Temperature, y = Q97.5, col = Species), 
            linetype = "dotted") +
  geom_line(data = d_Kelso, aes(x = Temperature, y = Estimate, col = Species), size = 1.5) +
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  geom_hline(yintercept = 0, col = "black") +
  scale_x_continuous( breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32),
                      labels = c(19, 21, 23, 25, "", 27, "", 29, "", 31, "") ) +
  xlab( "Temperature (°C)") +
  
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text.align = 0) +
  scale_y_continuous(breaks = seq(0, 6, by = 1)) +
  new_scale_fill() +
  geom_tile(data = intervals, aes(x = temp, y = y_val, fill = num_days), 
            height = 0.2) +
  scale_fill_viridis(option = "C", direction = -1,
                     name = "Days per year") 

Fig1


# Plot change in weight for the empty plugs compared to the estimated change in
# coral fragment weight

FigS5 <-
  ggplot() +
  
  geom_jitter(data = d[d$Reef == "Kelso", ], aes(x = Temperature, y = Growth, 
                                                 fill = Species),
              pch = 21, col = "black", size = 2, alpha = 0.3, height = 0,
              width = 0.2) +
  geom_line(data = d_Kelso, aes(x = Temperature, y = Estimate, col = Species), size = 1,
            alpha = 0.5) +
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  scale_x_continuous( breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32),
                      labels = c(19, 21, 23, 25, "", 27, "", 29, "", 31, "") ) +
  xlab( "Temperature (°C)") +
  geom_point(data = plugs, aes(x = Temperature, y = Growth),
             pch = 21, fill = "red", col = "black", size = 3, alpha = 0.5) +
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text.align = 0) +
  scale_y_continuous(breaks = seq(-1, 10, by = 1)) 
FigS5



# Plot raw data with model predictions
plot_raw <-
  ggplot() +
  geom_ribbon(data = d_Kelso,
              aes(x = Temperature, ymin = Q2.5, ymax = Q97.5, 
                  fill = Species), alpha = 0.2) +
  geom_line(data = d_Kelso,
            aes(x = Temperature, y = Q2.5, col = Species), 
            linetype = "dotted") +
  geom_line(data = d_Kelso,
            aes(x = Temperature, y = Q97.5, col = Species), 
            linetype = "dotted") +
  geom_line(data = d_Kelso,
            aes(x = Temperature, y = Estimate, col = Species), size = 1.5) +
  geom_jitter(data = d[d$Reef == "Kelso", ],
              aes(x = Temperature, y = Growth, fill = Species), pch = 21, 
              col ="black", width = 0.1, height = 0, alpha = 0.2) +
  facet_wrap( ~ Species) +
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  geom_hline(yintercept = 0, col = "black") +
  scale_x_continuous( breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32),
                      labels = c(19, 21, 23, 25, "", 27, "", 29, "", 31, "") ) +
  xlab( "Temperature (°C)") +
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text.align = 0,
        strip.text = element_text(face = "italic"), 
        strip.background = element_blank(),
        legend.position = "none")
plot_raw


Estimates <- as.data.frame(
  tidyMCMC(mod, conf.int = TRUE, conf.level = 0.95,
           conf.method = "HPDinterval")[1:21, ] )

#write.table(Estimates, "Estimates_TPC_model_Feb2022.csv", sep = ",", row.names = FALSE)



# Function to get thermal breadth and Topt
therm_br_fun <- function( rho, sigma, xshift,  prop) {
  Temps <- seq(0, 50, by = 0.01)
  
  prop_growth <- exp( - exp(rho * (Temps - xshift)) - sigma*(Temps - xshift)^ 2)/
    max(exp( - exp(rho * (Temps - xshift)) - sigma*(Temps - xshift)^ 2))
  Topt <- Temps[which(prop_growth == max(prop_growth))]
  
  n1 <- which(abs(prop_growth - prop) == min(abs(prop_growth - prop)))
  
  temp1 <- Temps[n1]
  
  Temps <- Temps[-n1]
  prop_growth <- prop_growth[-n1]
  
  n2 <- which(abs(prop_growth - prop) == min(abs(prop_growth - prop)))
  
  if ( temp1 <= Topt) {
    
    while (Temps[n2] <= Topt) {
      Temps <- Temps[-n2]
      prop_growth <- prop_growth[-n2]
      n2 <- which(abs(prop_growth - prop) == min(abs(prop_growth - prop)))
    }
    temp2 <- Temps[n2]
  }
  
  
  if ( temp1 >= Topt) {
    
    while (Temps[n2] >= Topt) {
      Temps <- Temps[-n2]
      prop_growth <- prop_growth[-n2]
      n2 <- which(abs(prop_growth - prop) == min(abs(prop_growth - prop)))
    }
    temp2 <- Temps[n2]
  }
  
  
  
  crit_temps <- c(temp1, temp2)
  
  return(list("Tbr_temps" = crit_temps[order(crit_temps)],
              "Topt" = Topt))
}





# Function to get the posterior distribution of thermal breadth and Topt
get_post_dist_Tbr <- function( Species, Reef, prop1, temps) {
  
  model <- mod
  
  
  if (Species == "Acropora hyacinthus") {
    post <- posterior_samples(model)[ , c(1, 5, 9, 13)]
  }
  if (Species == "Acropora tenuis") {
    post <- posterior_samples(model)[ , c(2, 6, 10, 14)]
  }   
  
  if (Species == "Pocillopora verrucosa") {
    post <- posterior_samples(model)[ , c(3, 7, 11, 15)]
  }
  
  if (Species == "Stylophora pistillata") {
    post <- posterior_samples(model)[ , c(4, 8, 12, 16)]
  }        
  
  
  
  
  df <- data.frame("T_min" = NA, "T_max" = NA, "num_of_days" = NA,
                   "Species" = NA,"Reef" = NA, "Gmax" = NA, "Topt" = NA)
  
  
  Gmax       <-      post[ ,1]
  rho        <-      post[ ,2]
  xshift     <-      post[ ,3]
  sigma      <-  exp(post[ ,4])
  
  T_min <- c()
  T_max <- c()
  Topt  <- c()
  
  num_of_days <- c()
  
  for (j in 1:length(Gmax)) {
    outputs <- therm_br_fun(rho        = rho[j], 
                            sigma      = sigma[j],
                            xshift     = xshift[j], 
                            prop       = prop1)
    
    T_min[j] <- outputs$Tbr_temps[1]
    T_max[j] <- outputs$Tbr_temps[2]
    Topt[j]     <- outputs$Topt
    
    if (is.na(T_min[j]) == FALSE & is.na(T_max[j]) == FALSE) {
      
      n <- which(temps >= T_min[j])
      temps1 <- temps[n]
      n <- which(temps1 <= T_max[j])
      temps1 <- temps1[n]
      num_of_days[j] <- (length(temps1) / length(temps)) * 365
      
    }
    
    
    
    
  }
  
  df2 <- data.frame(
    "T_min"       = T_min,
    "T_max"       = T_max,
    "num_of_days" = num_of_days,
    "Species"     = rep(Species,length(Gmax)),
    "Reef"        = rep(Reef, length(Gmax)),
    "Gmax"        = Gmax,
    "Topt"        = Topt)
  
  df <- rbind(df, df2)
  
  
  
  df$Tbr <- df$T_max - df$T_min
  
  df <- df[-1, ]
  
  
  sum_df <- df %>% 
    dplyr::group_by(Species, Reef) %>%
    dplyr::summarise(
      "Num_days"      = quantile(num_of_days, 0.500, na.rm = TRUE),
      "Num_days_lwr"  = quantile(num_of_days, 0.025, na.rm = TRUE),
      "Num_days_upr"  = quantile(num_of_days, 0.975, na.rm = TRUE),
      "Tbr"           = quantile(Tbr,         0.500, na.rm = TRUE),
      "Tbr_lwr"       = quantile(Tbr,         0.025, na.rm = TRUE),
      "Tbr_upr"       = quantile(Tbr,         0.975, na.rm = TRUE),
      "Tmin"          = quantile(T_min,       0.500, na.rm = TRUE),
      "Tmin_lwr"      = quantile(T_min,       0.025, na.rm = TRUE),
      "Tmin_upr"      = quantile(T_min,       0.975, na.rm = TRUE),
      "Tmax"          = quantile(T_max,       0.500, na.rm = TRUE),
      "Tmax_lwr"      = quantile(T_max,       0.025, na.rm = TRUE),
      "Tmax_upr"      = quantile(T_max,       0.975, na.rm = TRUE),
      "Gmax"          = quantile(Gmax,        0.500, na.rm = TRUE),
      "Gmax_lwr"      = quantile(Gmax,        0.025, na.rm = TRUE),
      "Gmax_upr"      = quantile(Gmax,        0.975, na.rm = TRUE),
      "Topt"          = quantile(Topt,        0.500, na.rm = TRUE),
      "Topt_lwr"      = quantile(Topt,        0.025, na.rm = TRUE),
      "Topt_upr"      = quantile(Topt,        0.975, na.rm = TRUE)
    )
  sum_df$Species <- Species
  sum_df$Reef <- Reef
  
  
  return(list("sum_df" = as.data.frame(sum_df), "all_values" = df))
  
  
}

# Get thermal breadth and Topt for each species
tbr_AH_K <- get_post_dist_Tbr("Acropora hyacinthus", "Kelso", 0.8, 
                              sub_temp$depth8m)
tbr_AT_K <- get_post_dist_Tbr("Acropora tenuis", "Kelso", 0.8, 
                              sub_temp$depth8m)
tbr_PV_K <- get_post_dist_Tbr("Pocillopora verrucosa", "Kelso", 0.8, 
                              sub_temp$depth8m)
tbr_SP_K <- get_post_dist_Tbr("Stylophora pistillata", "Kelso", 0.8, 
                              sub_temp$depth8m)


days_Species <- rbind(tbr_AH_K$sum_df, tbr_AT_K$sum_df, tbr_PV_K$sum_df,
                      tbr_SP_K$sum_df)
days_Species$y_val <- c(0, 0.2, 0.4, 0.6)



Fig1 <- Fig1 +
  geom_linerange(data = days_Species, aes(y = -0.4 - y_val, xmin = Tmin_lwr, xmax = Tmax_upr, 
                                          col = Species, group = Species), alpha = 0.8, size = 1,
                 show.legend = FALSE) +
  geom_linerange(data = days_Species, aes(y = -0.4 - y_val, xmin = Tmin, xmax = Tmax, 
                                          col = Species, group = Species), size = 3,
                 show_legend = FALSE) +
  geom_text(data = days_Species, aes(y = -0.4 - y_val, x = Tmax_upr + 0.6, 
                                     col = Species,
                                     label = round(Num_days)),
            show_legend = FALSE)
Fig1
#saveRDS(Fig1, "Fig1.rds")


# Get correlation between Gmax and Topt
all_pars <- c(
  "b_Gmax_SpeciesAcroporahyacinthus", "b_Gmax_SpeciesAcroporatenuis",
  "b_Gmax_SpeciesPocilloporaverrucosa", "b_Gmax_SpeciesStylophorapistillata"
)

n <- which(colnames(posterior_samples(mod)) %in% all_pars)
all_draws <- posterior_samples(mod)[ , n]
Topts <- cbind(tbr_AH_K$all_values$Topt, tbr_AT_K$all_values$Topt,
               tbr_PV_K$all_values$Topt, tbr_SP_K$all_values$Topt)
colnames(Topts) <- c("Topt_Ahyacinthus", "Topt_Atenuis", "Topt_Pverrucosa",
                     "Topt_Spistillata")
all_draws <- cbind(all_draws, Topts)

cor_df <- data.frame(r = NA, p = NA, draw = NA, intercept = NA, slope = NA)

for (j in 1:nrow(all_draws)) {
  df_cor <- data.frame("Gmax" = c(all_draws$b_Gmax_SpeciesAcroporahyacinthus[j],
                                  all_draws$b_Gmax_SpeciesAcroporatenuis[j],
                                  all_draws$b_Gmax_SpeciesPocilloporaverrucosa[j],
                                  all_draws$b_Gmax_SpeciesStylophorapistillata[j]),
                       "Topt" = c(all_draws$Topt_Ahyacinthus[j],
                                  all_draws$Topt_Atenuis[j],
                                  all_draws$Topt_Pverrucosa[j],
                                  all_draws$Topt_Spistillata[j]))
  m <- lm(df_cor$Gmax ~ df_cor$Topt)
  cor_j <- cor.test(df_cor$Gmax, df_cor$Topt)
  cor_df[j, ] <- c( cor_j$estimate,  cor_j$p.value,  j, coef(m)[1], coef(m)[2])
  
}


df_sum <- data.frame("Gmax" = c(all_draws$b_Gmax_SpeciesAcroporahyacinthus,
                                all_draws$b_Gmax_SpeciesAcroporatenuis,
                                all_draws$b_Gmax_SpeciesPocilloporaverrucosa,
                                all_draws$b_Gmax_SpeciesStylophorapistillata),
                     "Topt" = c(all_draws$Topt_Ahyacinthus,
                                all_draws$Topt_Atenuis,
                                all_draws$Topt_Pverrucosa,
                                all_draws$Topt_Spistillata),
                     "Species" = c(rep("Acropora hyacinthus", nrow(all_draws)),
                                   rep("Acropora tenuis", nrow(all_draws)),
                                   rep("Pocillopora verrucosa", nrow(all_draws)),
                                   rep("Stylophora pistillata", nrow(all_draws))))

quantile(cor_df$r, c(0.025, 0.975))


# Get quantiles of Gmax and Topt
summary_df <- data.frame("Species" = c("Acropora hyacinthus",
                                       "Acropora tenuis",
                                       "Pocillopora verrucosa",
                                       "Stylophora pistillata"),
                         "Gmax" = c(
                           quantile(all_draws$b_Gmax_SpeciesAcroporahyacinthus, 0.5),
                           quantile(all_draws$b_Gmax_SpeciesAcroporatenuis, 0.5),
                           quantile(all_draws$b_Gmax_SpeciesPocilloporaverrucosa, 0.5),
                           quantile(all_draws$b_Gmax_SpeciesStylophorapistillata, 0.5)
                         ),
                         "Gmax_lwr" = c(
                           quantile(all_draws$b_Gmax_SpeciesAcroporahyacinthus, 0.025),
                           quantile(all_draws$b_Gmax_SpeciesAcroporatenuis, 0.025),
                           quantile(all_draws$b_Gmax_SpeciesPocilloporaverrucosa, 0.025),
                           quantile(all_draws$b_Gmax_SpeciesStylophorapistillata, 0.025)
                         ),
                         "Gmax_upr" = c(
                           quantile(all_draws$b_Gmax_SpeciesAcroporahyacinthus, 0.975),
                           quantile(all_draws$b_Gmax_SpeciesAcroporatenuis, 0.975),
                           quantile(all_draws$b_Gmax_SpeciesPocilloporaverrucosa, 0.975),
                           quantile(all_draws$b_Gmax_SpeciesStylophorapistillata, 0.975)
                         ),
                         "Topt" = c(
                           quantile(all_draws$Topt_Ahyacinthus, 0.5),
                           quantile(all_draws$Topt_Atenuis, 0.5),
                           quantile(all_draws$Topt_Pverrucosa, 0.5),
                           quantile(all_draws$Topt_Spistillata, 0.5)
                         ),
                         "Topt_lwr" = c(
                           quantile(all_draws$Topt_Ahyacinthus, 0.025),
                           quantile(all_draws$Topt_Atenuis, 0.025),
                           quantile(all_draws$Topt_Pverrucosa, 0.025),
                           quantile(all_draws$Topt_Spistillata, 0.025)
                         ),
                         "Topt_upr" = c(
                           quantile(all_draws$Topt_Ahyacinthus, 0.975),
                           quantile(all_draws$Topt_Atenuis, 0.975),
                           quantile(all_draws$Topt_Pverrucosa, 0.975),
                           quantile(all_draws$Topt_Spistillata, 0.975)
                         )
)

summary_df <- as.data.frame(summary_df)

# Plot distribution of correlation coefficients
corr_hist <- ggplot(data = cor_df) +
  geom_histogram(mapping = aes(x = r), fill = "grey", colour = "grey60") +
  labs(x = expression(paste("Pearson's correlation between ", G[max], " and ",
                            T[opt])),
       y = "Frequency") +
  theme_bw() +
  geom_vline(aes( xintercept = quantile(r, 0.5)), col = "red", alpha = 0.5, 
             size = 1.2) +
  geom_vline(aes( xintercept = quantile(r, 0.025)), col = "red", alpha = 0.5, 
             size = 1, linetype = "dashed") +
  geom_vline(aes( xintercept = quantile(r, 0.975)), col = "red", alpha = 0.5, 
             size = 1, linetype = "dashed") +
  theme(text = element_text(size = 15))
corr_hist



min_Topt <- min(summary_df$Topt_lwr)
max_Topt <- max(summary_df$Topt_upr)

cor_lines <- matrix(NA, ncol = nrow(cor_df), nrow = 100)
Temp <- seq(min_Topt, max_Topt, l = 100)


for (i in 1:nrow(cor_df)) {
  cor_lines[ ,i] <- cor_df$intercept[i] + cor_df$slope[i] * Temp
}

cor_lines <- data.frame("Temperature" = Temp,
                        "lwr90" = apply(cor_lines, 1, function(x) quantile(x, 0.05)),
                        "upr90" = apply(cor_lines, 1, function(x) quantile(x, 0.95)),
                        "lwr80" = apply(cor_lines, 1, function(x) quantile(x, 0.1)),
                        "upr80" = apply(cor_lines, 1, function(x) quantile(x, 0.9)),
                        "lwr70" = apply(cor_lines, 1, function(x) quantile(x, 0.15)),
                        "upr70" = apply(cor_lines, 1, function(x) quantile(x, 0.85)),
                        "lwr60" = apply(cor_lines, 1, function(x) quantile(x, 0.2)),
                        "upr60" = apply(cor_lines, 1, function(x) quantile(x, 0.8)),
                        "lwr50" = apply(cor_lines, 1, function(x) quantile(x, 0.25)),
                        "upr50" = apply(cor_lines, 1, function(x) quantile(x, 0.75)),
                        "lwr40" = apply(cor_lines, 1, function(x) quantile(x, 0.3)),
                        "upr40" = apply(cor_lines, 1, function(x) quantile(x, 0.7)),
                        "lwr30" = apply(cor_lines, 1, function(x) quantile(x, 0.35)),
                        "upr30" = apply(cor_lines, 1, function(x) quantile(x, 0.65)),
                        "lwr20" = apply(cor_lines, 1, function(x) quantile(x, 0.4)),
                        "upr20" = apply(cor_lines, 1, function(x) quantile(x, 0.6)),
                        "lwr10" = apply(cor_lines, 1, function(x) quantile(x, 0.45)),
                        "upr10" = apply(cor_lines, 1, function(x) quantile(x, 0.55))
)



plot_4A <- ggplot() +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr90, ymax = upr90), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr80, ymax = upr80), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr70, ymax = upr70), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr60, ymax = upr60), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr50, ymax = upr50), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr40, ymax = upr40), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr30, ymax = upr30), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr20, ymax = upr20), 
              alpha = 0.05) +
  geom_ribbon(data = cor_lines, aes( x = Temperature, ymin = lwr10, ymax = upr10), 
              alpha = 0.05) +
  geom_linerange(data = summary_df,
                 aes(x = Topt, ymin = Gmax_lwr, ymax = Gmax_upr, col = Species),
                 size = 2, alpha = 0.6) +
  geom_linerangeh(data = summary_df,
                  aes(xmin = Topt_lwr, xmax = Topt_upr, y = Gmax, col = Species),
                  size = 2, alpha = 0.6) +
  geom_point(data = summary_df,
             aes(x = Topt, y = Gmax, fill = Species), pch = 21, col = "black", size = 5) +
  theme_bw()+
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  ylab(expression(paste(G[max]," (mg ", g^-1, d^-1, ")"))) +
  xlab(expression(paste(T[opt], " (°C)"))) +
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text.align = 0)
plot_4A

plot_4b <- ggplot(data = cor_df) +
  geom_histogram(mapping = aes(x = r), fill = "grey", colour = "grey60") +
  labs(x = "Correlation",
       y = "Frequency") +
  theme_bw() +
  geom_vline(aes( xintercept = quantile(r, 0.5)), col = "red", alpha = 0.5, 
             size = 1.2) +
  geom_vline(aes( xintercept = quantile(r, 0.025)), col = "red", alpha = 0.5, 
             size = 1, linetype = "dashed") +
  geom_vline(aes( xintercept = quantile(r, 0.975)), col = "red", alpha = 0.5, 
             size = 1, linetype = "dashed") +
  theme(text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        #panel.background = element_rect(col = "black"),
        legend.text.align = 0)
plot_4b

plot_4 <- plot_4A + gg_inset(ggplot2::ggplotGrob(plot_4b),
                             data = data.frame(group = 1),
                             xmin = 28.5, xmax= 30.4, 
                             ymin = 5, ymax = 7.5)
plot_4



# Get intraspecific variability
var_df <- data.frame("Parameter" = factor(), "Species" = factor(), 
                     "Colony" = factor(), "Var" = numeric(),
                     "Rel_var" = numeric())

post_samples <- posterior_samples(mod)

d <- d[!d$Reef == "No Name Reef", ]
d$colony <- factor(d$colony)

colonies <- unique(d$colony)


rand_rows <- sample(1 : nrow(post_samples), 1000, replace = FALSE)

for (i in rand_rows) {
  
  for (j in 1:length(colonies)) {
    
    Species <- paste(unlist(strsplit(as.character(colonies[j]), " "))[1], 
                     unlist(strsplit(as.character(colonies[j]), " "))[2], 
                     sep = " ")
    
    string_f <- str_replace( str_replace(colonies[j], " ", "."), " ", ".")
    string_f <- paste0(string_f, ",")
    n <- which(grepl(string_f, colnames(post_samples), fixed = FALSE) == TRUE)
    n_Sp <- which(grepl(sub(" ", "", Species), colnames(post_samples), 
                        fixed = FALSE) == TRUE)
    Gmax <- post_samples[i, n_Sp[1]]
    rho <- post_samples[i, n_Sp[2]]
    xshift <- post_samples[i, n_Sp[3]]
    sigma <- post_samples[i, n_Sp[4]]
    breadth <- therm_br_fun(rho = rho, sigma = exp(sigma),
                            xshift = xshift, 
                            prop = 0.8)
    Tbr <- abs(breadth$Tbr_temps[1] - breadth$Tbr_temps[2])
    Thigh <- breadth$Tbr_temps[2]
    Topt <- breadth$Topt
    
    if (xshift > 0) {
      breadth_colony <- therm_br_fun(rho = rho + post_samples[i, n[2]],
                                     xshift = xshift + post_samples[i , n[3]],
                                     sigma = exp(sigma + post_samples[i, n[4]]),
                                     prop = 0.8)
      Tbr_colony <- abs(breadth_colony$Tbr_temps[1] - breadth_colony$Tbr_temps[2])
      Topt_colony <- breadth_colony$Topt
      Thigh_colony <- breadth_colony$Tbr_temps[2]
      
    } else { 
      breadth_colony <- NA
      Tbr_colony <- NA 
      Topt_colony <- NA
      Thigh_colony <- NA
    }
    
    
    
    var_df_t <- data.frame("Parameter" = c("Gmax", "Topt", 
                                           "Tbr", "Thigh"),
                           "Species" = rep(Species, 4),
                           "Colony" = rep(colonies[j], 4),
                           "Var" = c(post_samples[i, n[1]],
                                     Topt - Topt_colony,
                                     Tbr - Tbr_colony,
                                     Thigh - Thigh_colony),
                           "Rel_var" = c(post_samples[i, n[1]]/Gmax,
                                         (Topt -Topt_colony)/Topt,
                                         (Tbr - Tbr_colony) / Tbr,
                                         (Thigh - Thigh_colony) / Thigh)
    )
    
    var_df <- rbind(var_df, var_df_t)
  }
  
  
}




#write.table(var_df, "Parameter_variation_modified_Frazier.csv", sep = ",",
#            row.names = FALSE)

var_df$Parameter <- as.factor(var_df$Parameter)
var_df$Colony <- as.factor(var_df$Colony)

summary_var <- var_df %>% 
  dplyr::group_by(Species, Colony, Parameter) %>%
  dplyr::summarise("median_var" = quantile(Var, 0.5, na.rm = TRUE),
                   "median_rel" = quantile(Rel_var, 0.5, na.rm = TRUE),
                   "upr_var" = quantile(Var, 0.975, na.rm = TRUE),
                   "upr_rel" = quantile(Rel_var, 0.975, na.rm =TRUE),
                   "lwr_var" = quantile(Var, 0.025, na.rm = TRUE),
                   "lwr_rel" = quantile(Rel_var, 0.025, na.rm = TRUE))

summary_var <- as.data.frame(summary_var)



# Plot estimate variability
plot_2B <-
  ggplot(summary_var, aes(x = Parameter, y = median_var , ymin = lwr_var ,
                          ymax = upr_var, col = Species, fill = Species)) +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(position = position_dodge2(width = 0.7), pch = 21, col = "black",
             alpha = 0.7) + 
  geom_linerange(position = position_dodge2(width = 0.7), alpha = 0.5) +
  scale_x_discrete(limits = c("Gmax", "Topt", "Tbr"),
                   labels = c(expression(paste(G[max], " (mg ", g^-1, d^-1, ")"),
                                         paste(T[opt], " (°C)"),
                                         paste(T[br], " (°C)")))) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  theme(text = element_text(size = 15),
        legend.text.align = 0,
        strip.text = element_text(face = "italic"), 
        strip.background = element_blank()) +
  ylab("Colony effects (+/- 95 % C.I.)") 
plot_2B


col_pred <- d[ , c("Species", "colony")]
col_pred <- col_pred[!duplicated(col_pred), ]
colonies_df <- data.frame("Species" = NA, "colony" = NA, "Temperature" = NA)
temp_vals <- seq(19, 31, by = 0.1)

for ( i in 1:nrow(col_pred)) {
  df <- data.frame("Species" = rep(col_pred$Species[i], length(temp_vals)),
                   "colony" = rep(col_pred$colony[i], length(temp_vals)),
                   "Temperature" = temp_vals)
  if (i == 1) {
    colonies_df <- df
  } else {
    colonies_df <- rbind(colonies_df, df)
  }
  
}

colonies_df <- cbind(colonies_df, as.data.frame(fitted(mod, colonies_df)))

col_pred <- col_pred[order(col_pred$Species), ]
col_pred$shapes <- rep(c(0:2, 5:6, 15:18), 4)

shapes <- col_pred$shapes
names(shapes) <- col_pred$colony

# Plot individual curves for each colony
plot_2A <-
  ggplot() +
  geom_point(data = d, aes(x = Temperature, y = Growth, col = Species, 
                           fill = Species, shape = colony), alpha = 0.3,
             show.legend = FALSE) + 
  geom_line(data = colonies_df, aes(x = Temperature, y = Estimate,
                                    col = Species, group = colony, alpha = 1),
            show.legend = FALSE) +
  scale_shape_manual(values = shapes) +
  scale_colour_manual(values = colours) +
  facet_wrap(~Species, ncol = 2) +
  ylab(expression(paste("Growth (mg ", g^-1, d^-1, ")"))) +
  scale_x_continuous( breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32),
                      labels = c(19, 21, 23, 25, "", 27, "", 29, "", 31, "") ) +
  xlab( "Temperature (°C)") + 
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text.align = 0,
        strip.text = element_text(face = "italic"), 
        strip.background = element_blank(),
        legend.position = "none") 

plot_2A


# Plot relative variation
plot_2B_supplementary <-
  ggplot(summary_var, aes(x = Parameter, y = median_rel * 100, ymin = lwr_rel * 100,
                          ymax = upr_rel * 100, col = Species, fill = Species)) +
  geom_hline(yintercept = 0, col = "black") +
  geom_point(position = position_dodge2(width = 0.7), pch = 21, col = "black",
             alpha = 0.7) + 
  geom_linerange(position = position_dodge2(width = 0.7), alpha = 0.5) +
  scale_x_discrete(limits = c("Gmax", "Topt", "Tbr"),
                   labels = c(expression(G[max], T[opt],
                                         T[br]))) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  theme_bw() +
  scale_colour_manual(values = colours, name = "", 
                      breaks = c("Acropora hyacinthus", "Acropora tenuis",
                                 "Pocillopora verrucosa", "Stylophora pistillata"),
                      labels = c(expression(italic("Acropora hyacinthus")),
                                 expression(italic("Acropora tenuis")),
                                 expression(italic("Pocillopora verrucosa")),
                                 expression(italic("Stylophora pistillata")))) +
  scale_fill_manual(values = colours, name = "", 
                    breaks = c("Acropora hyacinthus", "Acropora tenuis",
                               "Pocillopora verrucosa", "Stylophora pistillata"),
                    labels = c(expression(italic("Acropora hyacinthus")),
                               expression(italic("Acropora tenuis")),
                               expression(italic("Pocillopora verrucosa")),
                               expression(italic("Stylophora pistillata")))) +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(face = "italic"), 
        strip.background = element_blank(),
        legend.position = "none") +
  ylab("Percentage difference (+/- 95 % C.I.)") +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) 
plot_2B_supplementary 



Fig2 <- ggarrange(plot_2A + ggtitle('a'), plot_2B + ggtitle('b'))
Fig2


# Get table with estimates for Topt and Tbr
Estimates_derived <-
  data.frame("Variable" = c(rep("Topt", 4),
                           rep("Tbr", 4)),
             "Species" = rep(c("A. hyacinthus",
                               "A. tenuis",
                               "P. verrucosa",
                               "S. pistillata"), 2),
             "Estimate" = c(median(tbr_AH_K$all_values$Topt),
                            median(tbr_AT_K$all_values$Topt),
                            median(tbr_PV_K$all_values$Topt),
                            median(tbr_SP_K$all_values$Topt),
                            median(tbr_AH_K$all_values$Tbr),
                            median(tbr_AT_K$all_values$Tbr),
                            median(tbr_PV_K$all_values$Tbr),
                            median(tbr_SP_K$all_values$Tbr)),
             "S.D." = c(sd(tbr_AH_K$all_values$Topt),
                        sd(tbr_AT_K$all_values$Topt),
                        sd(tbr_PV_K$all_values$Topt),
                        sd(tbr_SP_K$all_values$Topt),
                        sd(tbr_AH_K$all_values$Tbr),
                        sd(tbr_AT_K$all_values$Tbr),
                        sd(tbr_PV_K$all_values$Tbr),
                        sd(tbr_SP_K$all_values$Tbr)),
             "Lower"= c(quantile(tbr_AH_K$all_values$Topt, 0.025),
                        quantile(tbr_AT_K$all_values$Topt, 0.025),
                        quantile(tbr_PV_K$all_values$Topt, 0.025),
                        quantile(tbr_SP_K$all_values$Topt, 0.025),
                        quantile(tbr_AH_K$all_values$Tbr, 0.025),
                        quantile(tbr_AT_K$all_values$Tbr, 0.025),
                        quantile(tbr_PV_K$all_values$Tbr, 0.025),
                        quantile(tbr_SP_K$all_values$Tbr, 0.025)),
             "Upper"= c(quantile(tbr_AH_K$all_values$Topt, 0.975),
                        quantile(tbr_AT_K$all_values$Topt, 0.975),
                        quantile(tbr_PV_K$all_values$Topt, 0.975),
                        quantile(tbr_SP_K$all_values$Topt, 0.975),
                        quantile(tbr_AH_K$all_values$Tbr, 0.975),
                        quantile(tbr_AT_K$all_values$Tbr, 0.975),
                        quantile(tbr_PV_K$all_values$Tbr, 0.975),
                        quantile(tbr_SP_K$all_values$Tbr, 0.975)))

#write.table(Estimates_derived, "Topt_Tbr_modified_Frazier.csv", sep = ",",
#            row.names = FALSE)
             