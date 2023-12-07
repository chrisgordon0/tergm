library(ggplot2)
library(threejs)
library(plotly)
library(patchwork)
library(ggiraph)
library(crosstalk)
library(dplyr)
library(tidyr)

#optpath <- read.csv('optpath.csv')
#optpath <- read.csv('/Users/chris/Documents/uni_work/tergm_data/EI-bayesOpt-200-10-200-50d-50iter.csv')
optpath <- read.csv('2dObjFunction.csv')
data <- data.frame(x=optpath$theta1, y=optpath$theta2, z=optpath$y)
#filtered_data <- data[data$z <= 1000, ]
filtered_data <- data
plot_ly(filtered_data, z=~log(1+z), y=~y, x=~x, type = "contour")

scatterplot3js(M, size=0.1, color=col, bg="black", pch=".")


optpath <- read.csv('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/log/EI-3.csv')
ggplot(optpath, aes(x=seq(1,length(y)),y=y)) + 
  geom_line() + 
  xlab("iteration") + 
  ylab("Mahalanobis Distance") +
  ggtitle("Mahalanobis Distance of target statistic over time")
  +


distances <- numeric()
for (i in 1:nrow(optpath)) {
  distances <- c(distances, sqrt((optpath[i,'theta1'] - (-5.202))^2 + (optpath[i,'theta2']-2.211)^2))
}
distances <- data.frame(dist=distances)
print(distances)

ggplot(distances, aes(x=seq(1,nrow(distances)), y=dist)) + 
  geom_line() + 
  xlab("iteration") + 
  ylab("Euclidian Distance") + 
  ggtitle("Convergence of Bayesian Optimization") + 
  geom_line(data=optpath, mapping=aes(x=seq(1,length(y)),y=y), colour='blue', alpha=0.6) + 
  theme_minimal()


################################################################################################################

plots <- lapply(1:5, function(i) {
  optpath <- read.csv(paste0('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOpt/EI/25/EI-', i, '.csv'))
  data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
  data[, 'des'] <- c(rep(1, 25), rep(0, 25))
  alpha <- c(rep(0.9, 25), seq(0.8, 1, length.out = 25))
  
  min_z_index <- which.min(data$z)
  data[,'min'] <- rep(0,50)
  data[min_z_index, 'min'] <- 1
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(colour = factor(des), alpha = alpha)) +
    theme_minimal() +
    geom_text(aes(label = ifelse(min,round(z,1),''), vjust = -0.6)) +
    scale_color_discrete(labels = c('Optimization', 'Design')) +
    labs(title = paste("Points Sampled During Bayesian Optimization using EI - File", i),
         colour = "Sampling Stage") +
    guides(alpha = "none", size = "none")
})

grid_plot <- wrap_plots(plots, ncol = 2)
grid_plot

### Label points near optimal

plots <- lapply(1:5, function(i) {
  optpath <- read.csv(paste0('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOpt/EI/25/EI-', i, '.csv'))
  data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
  data[, 'des'] <- c(rep(1, 25), rep(0, 25))
  alpha <- c(rep(0.9, 25), seq(0.8, 1, length.out = 25))
  
  min_z_index <- which.min(data$z)
  data[,'min'] <- rep(0, 50)
  data[min_z_index, 'min'] <- 1
  
  # Calculate Euclidean distance from (-5, 2)
  data$distance <- sqrt((data$x + 5)^2 + (data$y - 2)^2)
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(colour = factor(des), alpha = alpha)) +
    geom_text(data = data[data$distance <= 2 | data$min == 1, ], 
              aes(label = ifelse(min == 1 | (distance <= 2), round(z, 1), ''), vjust = -0.6), color = "black") +
    geom_point(data = data.frame(x = -5, y = 2), color = "black", size = 0.5) +  # Black dot at (-5, 2)
    theme_minimal() +
    scale_color_discrete(labels = c('Optimization', 'Design')) +
    labs(title = paste("Points Sampled During Bayesian Optimization using EI - File", i),
         colour = "Sampling Stage") +
    guides(alpha = "none", size = "none")
})

grid_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)
grid_plot




optpath <- read.csv('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOpt/EI/25/EI-1.csv')
data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)

plot_ly(data, z=~log(1+z), y=~y, x=~x, type = "contour")




######## AEI ########

plots <- lapply(1:5, function(i) {
  optpath <- read.csv(paste0('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/AEI/25/AEI-', i, '.csv'))
  data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
  data[, 'des'] <- c(rep(1, 25), rep(0, 25))
  alpha <- c(rep(0.9, 25), seq(0.8, 1, length.out = 25))
  
  min_z_index <- which.min(data$z)
  data[,'min'] <- rep(0,50)
  data[min_z_index, 'min'] <- 1
  
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(colour = factor(des), alpha = alpha)) +
    theme_minimal() +
    geom_text(aes(label = ifelse(min,round(z,1),''), vjust = -0.6)) +
    scale_color_discrete(labels = c('Optimization', 'Design')) +
    labs(title = paste("Points Sampled During Bayesian Optimization using AEI - File", i),
         colour = "Sampling Stage") +
    guides(alpha = "none", size = "none")
})

grid_plot <- wrap_plots(plots, ncol = 2)
grid_plot


optpath <- read.csv('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/AEI/25/AEI-3.csv')

data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
data[, 'des'] <- c(rep(1, 25), rep(0, 25))
alpha <- c(rep(0.9, 25), seq(0.8, 1, length.out = 25))
min_z_index <- which.min(data$z)
data[,'min'] <- rep(0, 50)
data[min_z_index, 'min'] <- 1

p <- ggplot(data, aes(x = x, y = y, text = sprintf("(%f,%f): %3f",x,y,z))) +
  geom_point(aes(colour = factor(des), alpha = alpha)) +
  theme_minimal() +
  geom_text(aes(label = ifelse(min, round(z, 1), ''), vjust = -0.6)) +
  scale_color_discrete(labels = c('Optimization', 'Design')) +
  labs(title = paste("Points Sampled During Bayesian Optimization using AEI - File",1),
       colour = "Sampling Stage") +
  guides(alpha = "none", size = "none")

ggplotly(p, tooltip = "text")


#######################

optpath <- read.csv('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/log/EI-1.csv')

data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
data[, 'des'] <- c(rep(1, 25), rep(0, 25))
min_z_index <- which.min(data$z)
data[,'min'] <- rep(0, 50)
data[min_z_index, 'min'] <- 1

p <- ggplot(data, aes(x = x, y = y, text = sprintf("(%f,%f): %3f",x,y,z))) +
  geom_point(aes(colour = factor(des))) +
  theme_minimal() +
  geom_text(aes(label = ifelse(min, round(z, 1), ''), vjust = -0.6)) +
  scale_color_discrete(labels = c('Optimization', 'Design')) +
  labs(title = paste("Points Sampled During Bayesian Optimization using EI - File",1),
       colour = "Sampling Stage") +
  guides(alpha = "none", size = "none")

ggplotly(p, tooltip = "text")





###### Average Convergence
distances <- numeric()
for (i in 1:30) {
  optpath <- read.csv(paste0('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/PowExp/log/UCB-PowExp-', i, '.csv'))
  data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
  min_z_index <- which.min(data$z)  
  distances <- c(distances, sqrt((optpath[min_z_index,'theta1'] - (-5.202))^2 + (optpath[min_z_index,'theta2']-2.211)^2))
  #print(optpath[min_z_index,'theta1'])
}
print(sum(distances<=0.1))
#distances
mean(distances)
var(distances)


############ Boxplot
distances <- list()

files <- c('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/EI-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/log/EI-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/AEI/25/AEI-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/AEI/25/log/AEI-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/UCB/25/UCB-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/UCB/25/log/UCB-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/CL/6x4/CL-',
           '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/CL/6x4/log/CL-')
names <- c('EI', 'EI+log', 'AEI', 'AEI+log', 'UCB', 'UCB+log', 'CL', 'CL+log')


# files <- c('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/log/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/sqrt/EI-')
# names <- c('EI', 'EI+log', 'EI+sqrt')

# files <- c('/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/exp+log/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/gauss+log/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/matern3_2+log/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/powexp+log/EI-',
#            '/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/EI/25/log/EI-')
# names <- c('Exponential', 'Gauss', 'Matérn(3,2)', 'Power Exponential', 'Matérn(5,2)')


for (fi in 1:length(files)) {
  dists <- numeric()
  for (i in 1:5) {
    optpath <- read.csv(paste0(files[fi], i, '.csv'))
    data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
    min_z_index <- which.min(data$z)  
    dists <- c(dists, sqrt((optpath[min_z_index,'theta1'] - (-5.202))^2 + (optpath[min_z_index,'theta2']-2.211)^2))
  }
  distances[[names[fi]]] <- dists
  
}

print(distances)
#mean(distances)
#var(distances)



distances <- data.frame(distances)

long_data <- distances %>%
  gather(key = "Category", value = "Value")

# Calculate the mean for each category
category_means <- long_data %>%
  group_by(Category) %>%
  summarise(mean_value = mean(Value))

# Sort the categories based on mean values
sorted_categories <- category_means %>%
  arrange(mean_value) %>%
  pull(Category)

# Reorder the Category factor levels based on sorted order
long_data$Category <- factor(long_data$Category, levels = sorted_categories)

# Create the box plot
plt <- ggplot(long_data, aes(x = Category, y = Value)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Method") +
  ylab("Euclidean distance to true parameter configuration")

plt

file_name <- "/Users/chris/Documents/uni_work/tergm_data/plots/EI_kernel_distances.png"

width <- 6
height <- 4

ggsave(file_name, plot = plt, width = width, height = height, dpi = 300)


##################### 2dBayesOptMore
############ Boxplot
distances <- list()


files <- c('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/plain/EI-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/AEI/plain/AEI-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/plain/UCB-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/plain/CL-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/log/EI-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/AEI/log/AEI-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/log/UCB-',
           '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/log/CL-')
names <- c('EI', 'AEI', 'UCB', 'CL', 'EI+log', 'AEI+log', 'UCB+log', 'CL+log')

# files <- c('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/Exp/log/EI-Exp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/Gauss/log/EI-Gauss-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/Mat32/log/EI-Mat32-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/PowExp/log/EI-PowExp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/EI/log/EI-')

# files <- c('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/Exp/log/CL-Exp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/Gauss/log/CL-Gauss-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/Mat32/log/CL-Mat32-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/PowExp/log/CL-PowExp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/CL/log/CL-')

# files <- c('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/Exp/log/UCB-Exp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/Gauss/log/UCB-Gauss-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/Mat32/log/UCB-Mat32-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/PowExp/log/UCB-PowExp-',
#            '/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/log/UCB-')
# 
# names <- c('Exponential', 'Gauss', 'Matérn(3,2)', 'Power Exponential', 'Matérn(5,2)')


for (fi in 1:length(files)) {
  dists <- numeric()
  for (i in 1:30) {
    optpath <- read.csv(paste0(files[fi], i, '.csv'))
    data <- data.frame(x = optpath$theta1, y = optpath$theta2, z = optpath$y)
    min_z_index <- which.min(data$z)  
    dists <- c(dists, sqrt((optpath[min_z_index,'theta1'] - (-5.202))^2 + (optpath[min_z_index,'theta2']-2.211)^2))
  }
  distances[[names[fi]]] <- dists
  
}

print(distances)
mean(distances[['EI']])
mean(distances[['AEI']])
mean(distances[['CL']])
mean(distances[['UCB']])
# 
# var(distances[['CL']])
# var(distances[['EI']])
# 

distances <- data.frame(distances)

long_data <- distances %>%
  gather(key = "Category", value = "Value")

# Calculate the mean for each category
category_means <- long_data %>%
  group_by(Category) %>%
  summarise(mean_value = mean(Value))

# Sort the categories based on mean values
sorted_categories <- category_means %>%
  arrange(mean_value) %>%
  pull(Category)

# Reorder the Category factor levels based on sorted order
long_data$Category <- factor(long_data$Category, levels = sorted_categories)

# Create the box plot
plt <- ggplot(long_data, aes(x = Category, y = Value)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Method") +
  ylab("Euclidean distance to true parameter configuration")

plt

file_name <- "/Users/chris/Documents/uni_work/tergm_data/new_plots/UCB_kernel_distances.png"

width <- 6
height <- 4

ggsave(file_name, plot = plt, width = width, height = height, dpi = 300)




#################### TIMES:
times <- read.csv('/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/times/plain/times.csv')#,  row.names=1)
times <- as.matrix(times, ncol=1, byrow=FALSE)
times
new_times <- c(times[30,1], times[,2], times[,3], times[,4])
new_times
# Convert hours to minutes (multiply by 60)
new_times <- c(new_times[1:32], new_times[33:length(new_times)] * 60)
new_times
# Calculate time intervals in minutes
time_intervals <- c(new_times[1],diff(new_times))
time_intervals

times <- c(times[1:29,1], time_intervals)
times <- matrix(times, ncol=4)

colnames(times) <- c('Gauss', 'Mat32', 'Exp', 'PowExp')
times
colMeans(times)
#write.csv(times, file='/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/times/times.csv',  row.names=FALSE)


times <- c(2.510632, 2.305852, 2.324209, 3.084200, 3.164384, 2.724171, 2.718787, 3.077633, 3.066180, 2.820555, 2.770765, 2.495589, 3.013746, 2.140674,2.314183, 2.159600, 2.364226, 2.180895, 2.510380, 2.844465, 2.597560, 2.966325, 2.604448, 2.830348, 2.709605, 2.969251, 2.717058, 2.457363, 2.580776, 2.435915, 4.997061, 7.836186, 10.270700, 13.384396, 16.147534, 18.460459, 21.129863, 23.422372, 25.924802, 28.081592, 31.680642, 33.909945, 36.573658, 39.350359, 41.959774, 44.625654, 47.646991, 50.721335, 53.317591, 56.043862, 58.912486, 1.026291, 1.074439, 1.111048, 1.158537, 1.202493, 1.244781, 1.281528, 1.324906, 1.372901, 1.411910, 1.463229, 1.507619, 1.555368, 1.602636, 1.648769, 1.690471, 1.735795, 1.784134, 1.831434, 1.873963, 1.909478, 1.954500, 2.003463, 2.047605, 2.087607, 2.139651, 2.184239, 2.225010, 2.274319, 2.321649, 2.365744, 2.410293, 2.453072, 2.500501, 2.544782, 2.583406, 2.621046, 2.667657, 2.709937, 2.758175, 2.806896, 2.853669, 2.897382, 2.942080, 2.996396, 3.043309, 3.086233, 3.131591, 3.174038, 3.220154, 3.266621, 3.307685, 3.351300, 3.394712, 3.441142, 3.497745, 3.541866, 3.586652, 3.644713, 3.692019, 3.738630, 3.791165, 3.834840, 3.878663, 3.928161, 3.967172, 4.012794, 4.053003, 4.099634)
new_times <- times[30:length(times)]
new_times <- c(new_times[1:22], new_times[23:length(new_times)] * 60)
time_intervals <- c(new_times[1],diff(new_times))
time_intervals
times <- c(times[1:29], time_intervals)
times <- matrix(times, ncol=4)
times
colnames(times) <- c('Gauss', 'Mat32', 'Exp', 'PowExp')
times
colMeans(times)
#write.csv(times, file='/Users/chris/Documents/uni_work/tergm_data_from_other/2dBayesOptMore/UCB/times/times.csv',  row.names=FALSE)




############## customOpt

results <- read.csv('/Users/chris/Documents/uni_work/tergm_data/customOptResults/resultsFromGrid.csv')
results
results[,'euclidean'] <- sqrt((-5.202-results[,3])^2+(2.197-results[,4])^2)
results[,'boolConverged'] <- sqrt((-5.202-results[,3])^2+(2.197-results[,4])^2) <= 0.2
ggplot(data=results, aes(y=euclidean)) + geom_boxplot() + theme_minimal()

results
#results[,6] <- results[,5]<0.1
sum(results[,'boolConverged'])
results
ind <- results[,'euclidean'] < 1
ind

results$custom_bins <- cut(results$euclidean, 
                           breaks = c(0, 0.5, 1.5, 2.5, Inf), 
                           labels = c("0-0.5", "0.5-1.5", "1.5-2.5", "2.5+"),
                           include.lowest = TRUE)

# To Plot The Grid of Starting Points:
plt1 <- ggplot(data=results, mapping=aes(x=X, y=X.1))+
  geom_point()+theme_minimal()+
  labs(x = "theta1", y = "theta2")+
  geom_point(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1), color = "red")+
  geom_text(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1, label = "true theta"), vjust = -0.5, color = "red")

ggsave('/Users/chris/Documents/uni_work/tergm_data/customOptResults/initGrid.png', plot = plt1, width = 6, height = 4, dpi = 300)


# Use viridis color palette for the custom bins
library(viridis)

plt <- ggplot(data = results, mapping = aes(x = X, y = X.1, color = custom_bins)) +
  geom_point() +
  geom_point(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1), color = "black") +
  geom_text(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1, label = "true theta"), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(x = "theta1", y = "theta2", color = "Converged") +
  scale_color_viridis(discrete = TRUE)

plt
ggsave('/Users/chris/Documents/uni_work/tergm_data/customOptResults/gridPlot0p2ppppppp.png', plot = plt, width = 6, height = 4, dpi = 300)



results <- read.csv('/Users/chris/Documents/uni_work/tergm_data/customOptResults/results.csv')
results
results[,'euclidean'] <- sqrt((-5.19-results[,3])^2+(2.10-results[,4])^2)
results[,'boolConverged'] <- sqrt((-5.19-results[,3])^2+(2.10-results[,4])^2) <= 0.2
ggplot(data=results, aes(y=euclidean)) + geom_boxplot() + theme_minimal()

results
#results[,6] <- results[,5]<0.1

results

plt <- ggplot(data = results, mapping = aes(x = X, y = X.1, color = boolConverged)) +
  geom_point() +
  geom_point(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1), color = "black") +
  geom_text(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1, label = "true theta"), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(x = "theta1", y = "theta2", color="Converged") + 
  scale_color_gradient(low = "blue", high = "red")

results

plt <- ggplot(data = results, mapping = aes(X, X.1, z=euclidean)) + 
  geom_contour_filled() + 
  theme_minimal() + 
  labs(x = "theta1", y = "theta2")


ggsave('/Users/chris/Documents/uni_work/tergm_data/customOptResults/contour.png', plot = plt, width = 6, height = 4, dpi = 300)


########################################
thetas <- read.csv('/Users/chris/Documents/uni_work/tergm_data/results/thetas.csv')
thetas

plt1 <- ggplot(data=thetas, mapping=aes(x=theta1, y=theta2))+
  geom_point()+theme_minimal()+
  labs(x = "theta1", y = "theta2")+
  geom_point(data = data.frame(X = -5.1, X.1 = 2.1), aes(x = X, y = X.1), color = "red")
plt1
ggsave('/Users/chris/Documents/uni_work/tergm_data/results/initConfigsAlgoComp.png', plot = plt1, width = 6, height = 4, dpi = 300)

bo <- read.csv('/Users/chris/Documents/uni_work/tergm_data/results/bayesOpt.csv')
co <- read.csv('/Users/chris/Documents/uni_work/tergm_data/results/customOpt.csv')
gd <- read.csv('/Users/chris/Documents/uni_work/tergm_data/results/GD.csv')

gd[31,'times'] <- gd[31,'times']/60
gd
mean(gd[,'times'])
var(gd[,'times'])

var(bo[,'times'])

co[which(co[,'times']>10),'times'] <- co[which(co[,'times']>10),'times']/60
co
mean(co[,'times'])
var(co[,'times'])
co
thetas
bo
bo[,'euclidean'] <- sqrt((-5.187386-bo[,1])^2+(2.197225-bo[,2])^2)
bo[,'category'] <- rep('Bayesian Optimization', times=50)
co[,'euclidean'] <- sqrt((-5.187386-co[,1])^2+(2.197225-co[,2])^2)
co[,'category'] <- rep('Alternate BO', times=50)
gd[,'euclidean'] <- sqrt((-5.187386-gd[,1])^2+(2.197225-gd[,2])^2)
gd[,'category'] <- rep('Gradient Descent', times=50)

sum(co[,'times'])/50

gd[which(gd$euclidean == Inf),'euclidean'] = 15


dat <- rbind(bo[,c('euclidean', 'category')], co[,c('euclidean', 'category')], gd[,c('euclidean', 'category')])

gd
thetas



p <- ggplot(data=dat, aes(x=category, y=euclidean)) + 
  geom_boxplot() + 
  theme_minimal()

labels_to_add <- data.frame(
  category = c('Gradient Descent'),
  value = c(15),
  label = c('Infinity (optimization failed)')
)

plt <- p + geom_text(data = labels_to_add, aes(x = category, y = value, label = label), vjust = 1.5) + xlab("Optimization Approach") + ylab("Euclidean Distance to True Parameter")

ggsave('/Users/chris/Documents/uni_work/tergm_data/results/comparison.png', plot = plt, width = 6, height = 4, dpi = 300)


thetas <- read.csv('/Users/chris/Documents/uni_work/tergm_data/results/thetas.csv')
ggplot(thetas, aes(x=theta1, y=theta2)) + theme_minimal() + geom_point()


