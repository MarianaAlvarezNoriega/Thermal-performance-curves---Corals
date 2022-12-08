

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

conv_dat <- read.csv("dry-wet_weight_conversion_coefficients.csv", sep = ",")

weeks <- read.csv("Sampling_weeks.csv", sep = ",")
colnames(weeks) <- c("Week", "Date")

weeks$Date <- as.Date(weeks$Date, format = "%d/%m/%Y")


dat <- read.csv("Buoyant_weight_Kelso_Reef.csv", sep = ",")
colnames(dat)[1] <- "Plug_number"

n <- grep("Date", colnames(dat))

dat[ ,n] <- lapply(dat[ , n], as.Date, format = "%d/%m/%Y")

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

dat2$Alive <- 1
dat2$Alive[grep("Sampled", dat2$Notes, ignore.case = TRUE)] <- 0
dat2$Alive[grep("Not Sampled", dat2$Notes, ignore.case = TRUE)] <- 1
dat2[dat2$position == 0, ]$Alive <- NA


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


dat2$Date <- as.factor(dat2$Date)
weeks$Date <- as.factor(weeks$Date)

dat2 <- merge(dat2, weeks, by = "Date", all.x = TRUE)


correct_week <- function(week, Temperature) {
  
  if (Temperature %in% c(19, 21, 23, 25, 31)) {
    
    w <- week - 2
  }
  
  else {
    w <- week - 1
  }
  return(w)
  
}

dat2$corrected_Week <- c()


for (i in 1:nrow(dat2)) {
  
  dat2$corrected_week[i] <- correct_week(dat2$Week[i], dat2$Temperature[i])
}




d <- data.frame("Species" = NA, "Reef" = NA, "Temperature" = NA,
                "Week" = NA, "Alive" = NA)


dat2 <- dat2[!dat2$Species == "plug", ]
Species <- unique(dat2$Species)
Temps <- unique(dat2$Temperature)

dat2 <- dat2[! is.na(dat2$Buoyant_weight) == T, ]

week_num <- seq(0, 8 , by = 2)

for (s in Species) {

  for (t in Temps) {
    
    for (w in week_num){
      sub_d <- dat2[dat2$Species == s & dat2$Reef == "Kelso", ]
      sub_d <- sub_d[sub_d$Temperature == t, ]
      sub_d <- sub_d[sub_d$corrected_week == w, ]
      
      surv <- data.frame("Species" = s,
                         "Reef" = "Kelso",
                         "Temperature" = t,
                         "Week" = w,
                         "Alive" = length(unique(sub_d$Plug))
                         
      )
      
      d <- rbind(d, surv)
      
    }
    
  }
}

d <- d[-1, ]

d$survival <- d$Alive / 9


p <-
  ggplot(d[d$Reef == "Kelso", ], aes(x = Week, y = Temperature, fill = survival)) +
  geom_tile(colour = "black") +
  theme_classic() +
  scale_fill_viridis(option = "C", direction = -1, name = "Proportion surviving") +
  facet_wrap(~Species) +
  scale_y_continuous(breaks = c(19, 21, 23, 25, 26, 27, 28, 29, 30, 31)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8)) +
  theme(strip.text = element_text(face = "italic")) +
  ylab( "Temperature (°C)") +
  xlab("Number of weeks") +
  theme(panel.background = element_rect(fill = "grey"),
        text = element_text(size = 15))
p  




