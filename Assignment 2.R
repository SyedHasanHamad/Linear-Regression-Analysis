library(faraway)
library(MASS)

data(fat)

n <- dim(fat)[1]

n1 <- round(n / 10)

flag <- c(1, 21, 22, 57, 70, 88, 91, 94, 121, 127, 149, 151, 159, 162,
          164, 177, 179, 194, 206, 214, 215, 221, 240, 241, 243)

fat1train <- fat[-flag, ]
fat1test <- fat[flag, ]


dim(fat1train)


hist(fat1train$brozek, main = "Histogram of Brozek's Equation", xlab = "Brozek's Equation")

summary(fat1train)

cor_matrix <- cor(fat1train[, -1])  
round(cor_matrix, 2)  





# Q (C) i. Linear Regression with all predictors
mod1 <- lm(brozek ~ ., data = fat1train)
# Linear Regression with all predictors
pred1.test <- predict.lm(mod1, newdata = fat1test)
error1.test <- mean((pred1.test - fat1test$brozek)^2)




library(leaps)
a2 <- regsubsets(brozek ~ ., data= fat1train, nvmax=5, nbest= 120,
method= c("exhaustive"), really.big= T)
models2 <- summary(a2)$which;
models2.size <- as.numeric(attr( models2, "dimnames")[[1]])
models2.rss <- summary(a2)$rss
op2 <- which( models2.size == 5)
flag2 <- op2[which.min( models2.rss[op2])]

mod2 <- lm( brozek ~ siri + density + thigh + knee +forearm, data = fat1train)
pred2 <- predict(mod2, fat1test[,2:18])
error2.test <- mean((pred2 - fat1test[,1])^2)


# Q (C) iii. Linear regression with variables (stepwise) selected using 
model3 <- step(mod1)
round(coef(model3),3)
pred3 <- predict(model3, fat1test[,2:18])
error3.test <- mean((pred3 - fat1test[,1])^2)


# q (C) iv. Ridge regression
library(MASS)
fat.ridge <- lm.ridge( brozek ~ ., data = fat1train, lambda= seq(0,100,0.001))
lambdaopt <- which.min(fat.ridge$GCV)
rig1coef <- fat.ridge$coef[,lambdaopt]
rig1intercepts <- fat.ridge$ym - sum(fat.ridge$xm *(rig1coef / fat.ridge$scales));
pred4 <- scale(fat1test[,2:18], center = F, scale = fat.ridge$scales)%*%
rig1coef + rig1intercepts;
error4.test <- mean( (pred4 - fat1test[,1])^2)


# q (C) v. LASSO regression
library(lars)
fat.lars <- lars( as.matrix(fat1train[,2:18]),
fat1train[,1], type= "lasso", trace= TRUE);
Cp1 <- summary(fat.lars)$Cp;
index1 <- which.min(Cp1);
lasso.lambda <- fat.lars$lambda[index1];
fit5b <- predict(fat.lars, fat1test[,-1], s=lasso.lambda, type="fit", mode="lambda");
yhat5b <- fit5b$fit;
error5.test <- mean( (yhat5b - fat1test[,1])^2);


trainpca <- prcomp(fat1train[, 2:18])
trainpca$sdev
round(trainpca$sdev, 2)
modelpca <- lm(brozek ~ trainpca$x[, 1:17], data = fat1train)

beta.pca <- trainpca$rot[, 1:17] %*% modelpca$coef[-1]
library(pls)
fat.pca <- pcr(brozek ~ ., data = fat1train, validation = "CV")

ncompopt <- which.min(fat.pca$validation$adj)
xmean <- apply(fat1train[, 2:18], 2, mean)
xtesttransform <- as.matrix(sweep(fat1test[, 2:18], 2, xmean))
xtestPC <- xtesttransform %*% trainpca$rot[, 1:17]
ypred6 <- cbind(1, xtestPC) %*% modelpca$coef

ypred6b <- predict(fat.pca, ncomp = ncompopt, newdata = fat1test[, 2:18])
error6.test <- mean((ypred6b - fat1test[, 1])^2)




library(pls)

fat.pls <- plsr(brozek ~ ., data = fat1train, validation = "CV")
ncompopt <- which.min(fat.pls$validation$adj)

ypred7b <- predict(fat.pls, fat1test[, 2:18], ncomp = ncompopt)
error7.test <- mean((ypred7b - fat1test[, 1])^2)




# Print testing errors rounded to 4 decimal places
cat("Testing Errors:\n")
print(paste("Linear Regression with all predictors:", round(error1.test, 4)))
print(paste("Linear regression with best subset k=5:", round(error2.test, 4)))
print(paste("Linear regression with variables selected using AIC:", round(error3.test, 5)))
print(paste("Ridge regression:", round(error4.test, 4)))
print(paste("LASSO regression:", round(error5.test, 4)))
print(paste("Principal Component regression:", round(error6.test, 4)))
print(paste("Partial least squares:", round(error7.test, 4)))




# Set the seed for reproducibility
set.seed(7406)

# Initialize lists to store testing errors for each model
mean_testing_errors_list <- list()
variance_testing_errors_list <- list()

# Monte Carlo Cross-Validation iterations
B <- 100

# Loop through each model
for (i in 1:7) {
  # Initialize vectors to store testing errors for MCCV
  mccv_errors <- numeric(B)
  
  for (b in 1:B) {
    # Randomly split the data into training and testing sets
    train_indices <- sample(1:nrow(fat1train), nrow(fat1train) * 0.9)  # Adjust the proportion as needed
    train_data <- fat1train[train_indices, ]
    test_data <- fat1train[-train_indices, ]
    
    if (i == 1) {
      # Linear Regression with all predictors
      mod <- lm(brozek ~ ., data = train_data)
      pred <- predict.lm(mod, newdata = test_data)
    } else if (i == 2) {
      # Linear Regression with 5 best predictors
      mod <- lm(brozek ~ siri + density + thigh + knee + forearm, data = train_data)
      pred <- predict(mod, newdata = test_data[, 2:18])
    } else if (i == 3) {
      # Linear regression with variables (stepwise) selected
      mod <- step(lm(brozek ~ ., data = train_data))
      pred <- predict(mod, newdata = test_data[, 2:18])
    } else if (i == 4) {
      # Ridge regression
      library(MASS)
      fat.ridge <- lm.ridge(brozek ~ ., data = train_data, lambda = seq(0, 100, 0.001))
      lambdaopt <- which.min(fat.ridge$GCV)
      rig1coef <- fat.ridge$coef[, lambdaopt]
      rig1intercepts <- fat.ridge$ym - sum(fat.ridge$xm * (rig1coef / fat.ridge$scales))
      pred <- scale(test_data[, 2:18], center = FALSE, scale = fat.ridge$scales) %*%
        rig1coef + rig1intercepts
    } else if (i == 5) {
      # LASSO regression
      library(lars)
      fat.lars <- lars(as.matrix(train_data[, 2:18]), train_data[, 1], type = "lasso", trace = TRUE)
      Cp1 <- summary(fat.lars)$Cp
      index1 <- which.min(Cp1)
      lasso.lambda <- fat.lars$lambda[index1]
      fit5b <- predict(fat.lars, test_data[, -1], s = lasso.lambda, type = "fit", mode = "lambda")
      pred <- fit5b$fit
    } else if (i == 6) {
      # Principal Component Analysis (PCA)
      trainpca <- prcomp(train_data[, 2:18])
      modelpca <- lm(brozek ~ trainpca$x[, 1:17], data = train_data)
      beta.pca <- trainpca$rot[, 1:17] %*% modelpca$coef[-1]
      fat.pca <- pcr(brozek ~ ., data = train_data, validation = "CV")
      ncompopt <- which.min(fat.pca$validation$adj)
      xmean <- apply(train_data[, 2:18], 2, mean)
      xtesttransform <- as.matrix(sweep(test_data[, 2:18], 2, xmean))
      xtestPC <- xtesttransform %*% trainpca$rot[, 1:17]
      pred <- cbind(1, xtestPC) %*% modelpca$coef
    } else if (i == 7) {
      # Partial Least Squares (PLS)
      library(pls)
      fat.pls <- plsr(brozek ~ ., data = train_data, validation = "CV")
      ncompopt <- which.min(fat.pls$validation$adj)
      pred <- predict(fat.pls, test_data[, 2:18], ncomp = ncompopt)
    }
    
    # Calculate testing error for the current MCCV iteration
    mccv_errors[b] <- mean((pred - test_data$brozek)^2)
  }
  
  # Calculate and store the mean MCCV error for the current model
  mean_mccv_error <- mean(mccv_errors)
  variance_mccv_error <- var(mccv_errors)
  mean_testing_errors_list[[i]] <- mean_mccv_error
  variance_testing_errors_list[[i]] <- variance_mccv_error
}




# Print the mean and variance for each model
for (i in 1:7) {
  cat("Model", i, ": Mean MCCV Error =", mean_testing_errors_list[[i]], ", Variance MCCV Error =", variance_testing_errors_list[[i]], "\n")
}



