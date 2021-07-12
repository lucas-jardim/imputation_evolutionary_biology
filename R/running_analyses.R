library("rpart")
library("rpart.plot")


# Run regression tree -----------------------------------------------------

k_data <- cbind(log(abs((matriz.1[, 2] - matriz.1[, 1])/matriz.1[, 1])), matriz.1[, 11:ncol(matriz.1)])
colnames(k_data)[1] <- "k_error"
k <- rpart(k_error ~ ., data = k_data)

moran_data <- cbind(log(abs((matriz.1[, 4] - matriz.1[, 3])/matriz.1[, 3])), matriz.1[, 11:ncol(matriz.1)])
colnames(moran_data)[1] <- "moran_error"
moran <- rpart(moran_error ~ ., data = moran_data)

beta_data <- cbind(log(abs((matriz.1[, 10] - matriz.1[, 9])/matriz.1[, 9])), matriz.1[, 11:ncol(matriz.1)])
colnames(beta_data)[1] <- "beta_error"
beta <- rpart(beta_error ~ ., data = beta_data)

mean_data <- cbind(log(abs((matriz.1[, 6] - matriz.1[, 5])/matriz.1[, 5])), matriz.1[, 11:ncol(matriz.1)])
colnames(mean_data)[1] <- "mean_error"
mean <- rpart(mean_error ~ ., data = mean_data)

var_data <- cbind(log(abs((matriz.1[, 8] - matriz.1[, 7])/matriz.1[, 7])), matriz.1[, 11:ncol(matriz.1)])
colnames(var_data)[1] <- "var_error"
var <- rpart(var_error ~ ., data = var_data)



# Extracting variable's importance ----------------------------------------

importance_variable <- matrix(0, nrow = 5, ncol = 5)
colnames(importance_variable) <- colnames(matriz.1)[11:ncol(matriz.1)]
row.names(importance_variable) <- c("k", "moran", "beta", "mean", "var")
results <- list(k, moran, beta, mean, var)

for(i in 1:5){
 vars <- results[[i]]$variable.importance/sum(results[[i]]$variable.importance)
 importance_variable[i, names(vars)] <- vars
}



