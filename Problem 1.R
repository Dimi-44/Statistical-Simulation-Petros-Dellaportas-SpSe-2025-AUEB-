#Check the shape of the density
# Creating the Sequence
gfg =seq(0, 1, by = 0.001)
 
# Plotting the beta density
plot(gfg, dbeta(gfg, 1.6,2.1), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red")
# Find the maximum of the Beta density function
alpha <- 1.6
beta <- 2.1
beta_density <- function(x) {
  dbeta(x, alpha, beta)
}
result <- optimize(beta_density, interval = c(0, 1), maximum = TRUE)
result$maximum;result$objective


#################################################################
################  Use Uniform as an envelope  ###################
#################################################################

plot(gfg, dbeta(gfg, 1.6,2.1), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red",main="Uniform Envelope")
abline(h=1.48 ,col="black",lty=2)
abline(v=0,col="black",lty=2)   
abline(v=1,col="black",lty=2)   


# Define the parameters of the Beta distribution
alpha <- 1.6
beta <- 2.1

# Define the maximum value of the Beta density function
M <- 1.48

# Function to sample from Beta(1.6, 2.1) using acceptance-rejection method
sample_beta <- function(n) {
  samples <- numeric(n)
  accepted <- 0
  total_proposed <- 0

  while (accepted < n) {
    # Sample x from U(0, 1)
    x <- runif(1)

    # Sample u from U(0, 1)
    u <- runif(1)

    # Calculate the acceptance ratio
    acceptance_ratio <- dbeta(x, alpha, beta) / M

    # Accept or reject the sample
    if (u <= acceptance_ratio) {
      accepted <- accepted + 1
      samples[accepted] <- x
    }

    total_proposed <- total_proposed + 1
  }

  # Calculate the acceptance probability
  acceptance_probability <- accepted / total_proposed

  list(samples = samples, acceptance_probability = acceptance_probability)
}

# Number of samples to generate
n <- 10000

# Generate samples and calculate acceptance probability
result <- sample_beta(n)

# Print the acceptance probability
cat("Acceptance Probability:", result$acceptance_probability, "\n")

# Print the first few samples
cat("First few samples:", head(result$samples), "\n")

#################################################################
################  Use Triangle as an envelope  ##################
#################################################################
# Load necessary libraries
library(ggplot2)

# Parameters
alpha <- 1.6
beta_param <- 2.1
mode <- (alpha - 1) / (alpha + beta_param - 2)  # Mode of Beta distribution
c <- 1.1  # Scaling factor for the envelope

# Generate x values
x <- seq(0, 1, length.out = 200)

# Beta distribution PDF
beta_dist <- dbeta(x, alpha, beta_param)

# Triangular envelope PDF
triangular_envelope <- function(x, mode, c) {
  ifelse(x < 0 | x > 1, 0,
         ifelse(x <= mode, c * (2 * x / mode), c * (2 * (1 - x) / (1 - mode))))
}

# Compute envelope values
envelope <- sapply(x, triangular_envelope, mode = mode, c = c)

# Combine data for plotting
data <- data.frame(x = x, beta_dist = beta_dist, envelope = envelope)

# Plotting
ggplot(data, aes(x = x)) +
  geom_line(aes(y = beta_dist, color = "Beta(1.6, 2.1) Density"), linewidth = 1) +
  geom_line(aes(y = envelope, color = paste0("Triangular Envelope (c=", c, ")")), linewidth = 1) +
  labs(x = "x", y = "Density", 
       title = "Triangular Envelope",
       color = "Legend") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("red", "blue")) +
  ylim(0, 2.5)  # Adjust y-limit to fit both curves

# Acceptance Probability Calculation
acceptance_probability <- 1 / c
cat(sprintf("Acceptance Probability: %.4f (or %.1f%%)\n", acceptance_probability, 100 * acceptance_probability))
#We reach the target acceptance probability but we don't cover the whole target distribution!!


#################################################################
################  Use Trapezium as an envelope  #################
#################################################################
# Set seed for reproducibility
set.seed(123)

# Parameters for Beta distribution
alpha <- 1.6
beta_param <- 2.1

# Parameters for the trapezoidal envelope
a <- 0.05005
b <- 0.55
h <- 1.48

# Area under the trapezoidal envelope h(x)
area_trapezoid <- (0.5 * a * h) + ((b - a) * h) + (0.5 * (1 - b) * h)  # 1.1096
c <- area_trapezoid  # Scaling factor

# Normalized trapezoidal density g(x) = h(x) / area
g <- function(x, a, b, h, area_trapezoid) {
  if (0 <= x && x <= a) {
    return(((h / a) * x) / area_trapezoid)
  } else if (a < x && x <= b) {
    return(h / area_trapezoid)
  } else if (b < x && x <= 1) {
    return((h * (1 - x) / (1 - b)) / area_trapezoid)
  } else {
    return(0)
  }
}

# Trapezoidal envelope h(x) = c * g(x)
h_func <- function(x, a, b, h) {
  if (0 <= x && x <= a) {
    return((h / a) * x)
  } else if (a < x && x <= b) {
    return(h)
  } else if (b < x && x <= 1) {
    return(h * (1 - x) / (1 - b))
  } else {
    return(0)
  }
}

# CDF of the normalized trapezoidal density g(x) for sampling
G <- function(x, a, b, h, area_trapezoid) {
  if (x <= 0) return(0)
  if (x >= 1) return(1)
  if (0 < x && x <= a) {
    return(0.5 * (h / a) * (x^2) / area_trapezoid)
  } else if (a < x && x <= b) {
    return((0.5 * h * a + h * (x - a)) / area_trapezoid)
  } else {
    return((0.5 * h * a + h * (b - a) + h * (x - b) - 0.5 * (h / (1 - b)) * (x - b)^2) / area_trapezoid)
  }
}

# Inverse CDF for sampling from g(x)
inverse_G <- function(u, a, b, h, area_trapezoid) {
  if (u < 0 || u > 1) stop("u must be between 0 and 1")
  
  # Area under each segment of g(x)
  area1 <- G(a, a, b, h, area_trapezoid)  # Area from 0 to a
  area2 <- G(b, a, b, h, area_trapezoid) - G(a, a, b, h, area_trapezoid)  # Area from a to b
  area3 <- 1 - G(b, a, b, h, area_trapezoid)  # Area from b to 1
  
  # Normalize u to the correct segment
  if (u <= area1) {
    # First segment (0 to a)
    return(sqrt(2 * u * area_trapezoid * a / h))
  } else if (u <= area1 + area2) {
    # Second segment (a to b)
    u_adjusted <- (u - area1) / area2
    return(a + u_adjusted * (b - a))
  } else {
    # Third segment (b to 1)
    u_adjusted <- (u - (area1 + area2)) / area3
    # Quadratic equation: -0.5 * (h/(1-b)) * x^2 + h * x - (term1 + u_adjusted * area_trapezoid) = 0
    term1 <- 0.5 * h * a + h * (b - a)
    A <- -0.5 * (h / (1 - b))
    B <- h
    C <- -(term1 + u_adjusted * area_trapezoid)
    discriminant <- B^2 - 4 * A * C
    if (discriminant < 0) {
      warning("Negative discriminant encountered, returning NA")
      return(NA)
    }
    x <- (-B + sqrt(discriminant)) / (2 * A)
    if (x < b || x > 1) {
      warning("Invalid x value, adjusting to boundary")
      return(1)  # Fallback to upper bound
    }
    return(x)
  }
}

# Rejection sampling function
rejection_sampling <- function(n, a, b, h, area_trapezoid) {
  samples <- numeric(n)
  accepted <- 0
  total_attempts <- 0
  
  while (accepted < n) {
    # Step 1: Sample x from g(x) using inverse CDF
    u1 <- runif(1)
    x <- inverse_G(u1, a, b, h, area_trapezoid)
    if (is.na(x)) next  # Skip if x is NA due to invalid discriminant
    
    # Step 2: Compute the ratio f(x) / (c * g(x))
    f_x <- dbeta(x, alpha, beta_param)  # Beta density
    h_x <- h_func(x, a, b, h)  # Trapezoidal envelope
    if (h_x == 0) next  # Avoid division by zero
    ratio <- f_x / h_x  # Since h(x) = c * g(x), this is f(x) / (c * g(x))
    
    # Step 3: Accept or reject
    u2 <- runif(1)
    total_attempts <- total_attempts + 1
    if (u2 <= ratio) {
      accepted <- accepted + 1
      samples[accepted] <- x
    }
  }
  
  # Empirical acceptance rate
  emp_acceptance_rate <- n / total_attempts
  return(list(samples = samples, emp_acceptance_rate = emp_acceptance_rate))
}

# Generate samples
n_samples <- 10000
result <- rejection_sampling(n_samples, a, b, h, area_trapezoid)
emp_acceptance_rate <- result$emp_acceptance_rate 

# Theoretical acceptance probability
theoretical_acceptance <- 1 / area_trapezoid
cat("Theoretical Acceptance Probability:", theoretical_acceptance, "\n")
cat("Empirical Acceptance Rate:", mean(emp_acceptance_rate), "\n")

# Plot the results
# Histogram of samples with true Beta density overlaid
hist(samples, breaks = 50, freq = FALSE, col = "lightblue", 
     main = "Rejection Sampling from Beta(1.6, 2.1) using Trapezoidal Envelope",
     xlab = "x", ylab = "Density")
curve(dbeta(x, alpha, beta_param), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Sampled Distribution", "True Beta(1.6, 2.1)"),
       fill = c("lightblue", "red"))

# Plot the trapezoidal envelope and Beta density
x <- seq(0, 1, length.out = 200)
beta_density <- dbeta(x, alpha, beta_param)
envelope <- sapply(x, function(xi) h_func(xi, a, b, h))

plot(x, beta_density, type = "l", col = "red", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "Trapezoidal Envelope")
lines(x, envelope, col = "black", lty = 2, lwd = 2)
legend("topright", legend = c("Beta(1.6, 2.1) Density", "Trapezoidal Envelope"),
       col = c("red", "black"), lty = c(1, 2), lwd = 2)
#It reaches the target probability theoritically but empirical sampling has 0.88%.
#Hence we are close but not at the desired level

#################################################################
################  Use Polygonal as an envelope  #################
#################################################################
# Set seed for reproducibility
set.seed(123)

# Parameters for Beta distribution
alpha <- 1.6
beta_param <- 2.1

# Define the polygonal envelope points
x_points <- c(0, 0.005, 0.050, 0.100, 0.20, 0.35, 0.50, 0.6, 0.80, 1)
h_points <- c(0, 0.443, 0.705, 1.086, 1.35, 1.54, 1.4, 1.25, 0.72, 0)

# Compute slopes for each segment
slopes <- diff(h_points) / diff(x_points)

# Polygonal envelope h(x)
h_func <- function(x, x_points, h_points, slopes) {
  if (x < 0 || x > 1) return(0)
  if (x == 0) return(h_points[1])
  if (x == 1) return(h_points[length(h_points)])
  
  # Find the segment
  for (i in 1:(length(x_points) - 1)) {
    if (x_points[i] <= x && x <= x_points[i + 1]) {
      return(h_points[i] + slopes[i] * (x - x_points[i]))
    }
  }
  return(0)  # Should not reach here
}

# Compute the area under the polygonal envelope
compute_area <- function(x_points, h_points) {
  area <- 0
  for (i in 1:(length(x_points) - 1)) {
    # Area of trapezoid for each segment
    area <- area + 0.5 * (h_points[i] + h_points[i + 1]) * (x_points[i + 1] - x_points[i])
  }
  return(area)
}

area_trapezoid <- compute_area(x_points, h_points)  # Should be ~1.1111375
c <- area_trapezoid  # Scaling factor

# Normalized polygonal density g(x) = h(x) / area
g <- function(x, x_points, h_points, slopes, area_trapezoid) {
  return(h_func(x, x_points, h_points, slopes) / area_trapezoid)
}

# CDF of the normalized polygonal density g(x) for sampling
G <- function(x, x_points, h_points, slopes, area_trapezoid) {
  if (x <= 0) return(0)
  if (x >= 1) return(1)
  
  cumulative_area <- 0
  for (i in 1:(length(x_points) - 1)) {
    if (x < x_points[i + 1]) {
      # Partial area up to x in the current segment
      x_start <- max(x_points[i], 0)
      x_end <- min(x, x_points[i + 1])
      h_start <- h_func(x_start, x_points, h_points, slopes)
      h_end <- h_func(x_end, x_points, h_points, slopes)
      segment_area <- 0.5 * (h_start + h_end) * (x_end - x_start)
      return(cumulative_area + segment_area / area_trapezoid)
    } else {
      # Full area of the segment
      segment_area <- 0.5 * (h_points[i] + h_points[i + 1]) * (x_points[i + 1] - x_points[i])
      cumulative_area <- cumulative_area + segment_area / area_trapezoid
    }
  }
  return(cumulative_area)
}

# Inverse CDF for sampling from g(x)
inverse_G <- function(u, x_points, h_points, slopes, area_trapezoid) {
  if (u < 0 || u > 1) stop("u must be between 0 and 1")
  
  # Compute cumulative areas at each x_point
  cumulative_areas <- numeric(length(x_points))
  for (i in 1:length(x_points)) {
    cumulative_areas[i] <- G(x_points[i], x_points, h_points, slopes, area_trapezoid)
  }
  
  # Find the segment where u lies
  for (i in 1:(length(x_points) - 1)) {
    if (u <= cumulative_areas[i + 1]) {
      # Interpolate within the segment
      u_start <- cumulative_areas[i]
      u_end <- cumulative_areas[i + 1]
      x_start <- x_points[i]
      x_end <- x_points[i + 1]
      h_start <- h_points[i]
      h_end <- h_points[i + 1]
      # Area up to x: 0.5 * (h_start + h(x)) * (x - x_start) = (u - u_start) * area_trapezoid
      slope <- slopes[i]
      A <- 0.5 * slope
      B <- h_start
      C <- -(u - u_start) * area_trapezoid
      discriminant <- B^2 - 4 * A * C
      if (discriminant < 0) {
        warning("Negative discriminant, returning boundary value")
        return(x_end)
      }
      t <- (-B + sqrt(discriminant)) / (2 * A)
      x <- x_start + t
      if (x < x_start || x > x_end) {
        warning("x out of bounds, returning boundary value")
        return(x_end)
      }
      return(x)
    }
  }
  return(1)  # Fallback to upper bound
}

# Rejection sampling function
rejection_sampling <- function(n, x_points, h_points, slopes, area_trapezoid) {
  samples <- numeric(n)
  accepted <- 0
  total_attempts <- 0
  
  while (accepted < n) {
    # Step 1: Sample x from g(x) using inverse CDF
    u1 <- runif(1)
    x <- inverse_G(u1, x_points, h_points, slopes, area_trapezoid)
    if (is.na(x)) next  # Skip if x is NA
    
    # Step 2: Compute the ratio f(x) / (c * g(x))
    f_x <- dbeta(x, alpha, beta_param)  # Beta density
    h_x <- h_func(x, x_points, h_points, slopes)  # Polygonal envelope
    if (h_x == 0) next  # Avoid division by zero
    ratio <- f_x / h_x  # Since h(x) = c * g(x), this is f(x) / (c * g(x))
    
    # Step 3: Accept or reject
    u2 <- runif(1)
    total_attempts <- total_attempts + 1
    if (u2 <= ratio) {
      accepted <- accepted + 1
      samples[accepted] <- x
    }
  }
  
  # Empirical acceptance rate
  emp_acceptance_rate <- n / total_attempts
  return(list(samples = samples, emp_acceptance_rate = emp_acceptance_rate))
}

# Generate samples
n_samples <- 10000
result <- rejection_sampling(n_samples, x_points, h_points, slopes, area_trapezoid)
samples <- result$samples
emp_acceptance_rate <- result$emp_acceptance_rate

# Theoretical acceptance probability
theoretical_acceptance <- 1 / area_trapezoid
cat("Theoretical Acceptance Probability:", theoretical_acceptance, "\n")
cat("Empirical Acceptance Rate:", emp_acceptance_rate, "\n")

# Plot the results
# Histogram of samples with true Beta density overlaid
hist(samples, breaks = 50, freq = FALSE, col = "lightblue", 
     main = "Rejection Sampling from Beta(1.6, 2.1) using Polygonal Envelope",
     xlab = "x", ylab = "Density")
curve(dbeta(x, alpha, beta_param), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Sampled Distribution", "True Beta(1.6, 2.1)"),
       fill = c("lightblue", "red"))

# Plot the polygonal envelope and Beta density
x <- seq(0, 1, length.out = 200)
beta_density <- dbeta(x, alpha, beta_param)
envelope <- sapply(x, function(xi) h_func(xi, x_points, h_points, slopes))

plot(x, beta_density, type = "l", col = "red", lwd = 2,
     xlab = "x", ylab = "Density",
     main = "Polygonal Envelope")
lines(x, envelope, col = "black", lty = 2, lwd = 3)
legend("topright", legend = c("Beta(1.6, 2.1) Density", "Polygonal Envelope"),
       col = c("red", "black"), lty = c(1, 2), lwd = 2)
#The acceptance probability is 97.1%. We will keep that solution 
#The idea of using polygonal is based on Karlis,Xekalaki work for the book
#"Advances in Mathematical and Statistical Modeling" of Arnold et.al