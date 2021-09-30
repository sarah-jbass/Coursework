# code prep for hw6

# clear workspace
rm(list=ls())
setwd("C:/z_toshiba/course work/phd/econ 761/hw/hw6/")

library(readr)
bidder_data <- read_table2("fpa.dat", col_names = FALSE)
colnames(bidder_data) <- c("bidder_1", "bidder_2", "bidder_3", "bidder_4")

bidder_data2 <- stack(bidder_data)
colnames(bidder_data2) <- c("bid", "bidder")

# density for each bidder
b1_pdf <- density(bidder_data$bidder_1, kernel = "epanechnikov")
b2_pdf <- density(bidder_data$bidder_2, kernel = "epanechnikov")
b3_pdf <- density(bidder_data$bidder_3, kernel = "epanechnikov")
b4_pdf <- density(bidder_data$bidder_4, kernel = "epanechnikov")
#plot(b1_pdf)
#plot(b2_pdf)
#plot(b3_pdf)
#plot(b4_pdf)

# cdf for each bidder
b1_cdf <- ecdf(bidder_data$bidder_1)
b2_cdf <- ecdf(bidder_data$bidder_2)
b3_cdf <- ecdf(bidder_data$bidder_3)
b4_cdf <- ecdf(bidder_data$bidder_4)
#plot(b1_cdf)
#plot(b2_cdf)
#plot(b3_cdf)
#plot(b4_cdf)

# density for all bidders
all_b_pdf <- density(bidder_data2$bid, kernel = "epanechnikov")
#plot(all_b_pdf)

# cdf for all bidders
all_b_cdf <- ecdf(bidder_data2$bid)
#plot(all_b_cdf)

# estimated values of Fu
all_b_cdf(quantile(bidder_data$bidder_1)[2])
all_b_cdf(quantile(bidder_data$bidder_2)[2])
all_b_cdf(quantile(bidder_data$bidder_3)[2])
all_b_cdf(quantile(bidder_data$bidder_4)[2])

all_b_cdf(quantile(bidder_data$bidder_1)[4])
all_b_cdf(quantile(bidder_data$bidder_2)[4])
all_b_cdf(quantile(bidder_data$bidder_3)[4])
all_b_cdf(quantile(bidder_data$bidder_4)[4])

# correlation between values
val_ind <- lm(bidder_data$bidder_1 ~ bidder_data$bidder_2+
                  bidder_data$bidder_3+bidder_data$bidder_4)
coefs <- summary(val_ind)$coefficients
coefs
