## fitting a linear regression model on COPD data ## 

COPD <- read.csv('CO.csv')

# fitting model #
# MWT1Best = outcome variable, FEV1 = predictor variable # 

MWT1Best_FEV1 <- lm(MWT1Best~FEV1, data = COPD)

# summary of model #

summary(MWT1Best_FEV1)

# plots # 
# plot 1 = checking homogeneity of variance/linear relation. no pattern = assumptions met #
# plot 2 = Q-Q plot, checks that residuals follow  normal dist. should fall on line # 

plot(MWT1Best_FEV1)


