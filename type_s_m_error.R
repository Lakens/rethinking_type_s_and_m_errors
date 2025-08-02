##########################################

# Function from Gelman & Carlin
# I added the unconditional Type S to the output - DL

df <- Inf
n.sims <- 10000
A <- 0.1
s <- 1
retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
  return(list(power=power, typeS=typeS, typeSunconditional=p.lo,exaggeration=exaggeration))
}

# Example: true effect size of 0.1, standard error 3.28, alpha=0.05
retrodesign(.1, 3.28)


# Compute Type 3 or S error
ncp <- (0.2 * sqrt(40)/2) # For t-distribution
crit_t <- qt(0.05 / 2, (20 * 2) - 2, 0)
pt(crit_t, 38, ncp)

# We can check against G*Power. We compute power for 2-sided test at 5% alpha,
# n = 20 per group, d = 0.2. We get 0.0945673. For a one-sided test at 2.5% alpha,
# we get power of 0.0895807. The difference is the Type 3 error.
0.0945673 - 0.0895807

# d = 0.2, sd = 1, alpha 0.05

retrodesign(0, 1/sqrt(19), 0.05)



# Example: true effect size of 0.1, standard error 3.28, alpha=0.05
res <- retrodesign(0.3, 1)
res
res$typeSunconditional/res$power
(res$typeSunconditional+res$power)
res$typeSunconditional/(res$typeSunconditional+res$power)



fprp <- function(alpha, power, prior) {
  # Calculate false positive probability
  false_positive <- (1 - prior) * alpha

  # Calculate true positive probability
  true_positive <- prior * power

  # Compute FPRP
  fprp_value <- false_positive / (false_positive + true_positive)

  return(fprp_value)
}

# Example: alpha = 0.05, power = 0.80, prior = 0.10
fprp(alpha = 0.025, power = 0.05, prior = 0.05)

# Example: true effect size of 2, standard error 8.1, alpha=0.05
retrodesign(0.001, 1/sqrt(20))
# Producing Figures 2a and 2b for the Gelman and Carlin paper
D_range <- c(seq(0,1,.01),seq(1,10,.1),100)
n <- length(D_range)
power <- rep(NA, n)
typeS <- rep(NA, n)
exaggeration <- rep(NA, n)
for (i in 1:n){
  a <- retrodesign(D_range[i], 1)
  power[i] <- a$power
  typeS[i] <- a$typeS
  exaggeration[i] <- a$exaggeration
}




retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000){
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
  return(list(power=power, typeS=typeS, exaggeration=exaggeration))
}

# Example: true effect size of 0.1, standard error 3.28, alpha=0.05
retrodesign(5.3, 3.28)
# Example: true effect size of 2, standard error 8.1, alpha=0.05
retrodesign(2, 8.1)


# Pearson's correlation
PRDA::retrospective(effect_size = .4, sample_n1 = 13, test_method = "pearson",
                    alternative = c("two_sided"),
                    sig_level = 0.05, B = 1000000)

# Two-sample t-test
PRDA::retrospective(effect_size = .4, sample_n1 = 20, sample_n2 = 20,
              test_method = "two_sample")
# Welch t-test
PRDA::retrospective(effect_size = .4, sample_n1 = 20, sample_n2 = 20,
              test_method = "welch", ratio_sd = 1.5)
# Paired t-test
PRDA::retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = 25,
              test_method = "paired")
# One-sample t-test
retrospective(effect_size = .3, sample_n1 = 25, sample_n2 = NULL,
              test_method = "one_sample")


library(BUCSS)
ss.power.it(t.observed=2.4, n=200, alpha.prior=.05, alpha.planned=.05,
            assurance=.50, power=.80, step=.001)



# CHatGPT code to correct effect sizes using Taylor and Muller, 1996

# Install needed packages if not already installed
# install.packages(c("effsize", "weightr", "puniform"))

library(effsize)   # for Cohen's d
library(weightr)   # for selection model bias correction
library(puniform)  # for p-curve bias correction

# Step 1: Simulate two groups (independent samples)
set.seed(123)

n <- 20  # per group
mu1 <- 0
mu2 <- 0.7  # true effect size
sd <- 1

group1 <- rnorm(n, mean = mu1, sd = sd)
group2 <- rnorm(n, mean = mu2, sd = sd)

# Step 2: Compute Cohen's d from raw data
cohen_result <- cohen.d(group1, group2, pooled = TRUE)
observed_d <- as.numeric(cohen_result$estimate)
cat("Observed Cohen's d:", observed_d, "\n")

# Step 3: Compute SE for d assuming equal group sizes
compute_se_d <- function(n_per_group) {
  sqrt(2 / n_per_group)
}
se_d <- compute_se_d(n)

# Step 4: Apply Taylor & Muller correction
correct_effect_size <- function(observed_d, se, alpha = 0.05) {
  c <- qnorm(1 - alpha / 2)
  objective <- function(delta) {
    bias <- (dnorm(c - delta / se) - dnorm(-c - delta / se)) /
      (pnorm(delta / se - c) + pnorm(-delta / se - c)) * se
    (observed_d - bias - delta)^2
  }
  optimize(objective, interval = c(-2, 2))$minimum
}

corrected_tm <- correct_effect_size(observed_d, se_d)
cat("Taylor & Muller corrected d:", corrected_tm, "\n")

# Step 5: Use weightr for selection model correction
pval <- 2 * (1 - pnorm(abs(observed_d / se_d)))
dat_weightr <- data.frame(effect = observed_d, v = se_d^2, pval = pval)
dat_weightr$p_group <- cut(dat_weightr$pval, breaks = c(0, 0.025, 1), right = FALSE)

# Note: weightr requires multiple studies; this is for illustration only
weightr_result <- weightfunct(effect = dat_weightr$effect,
                              v = dat_weightr$v,
                              steps = c(0.025, 1),
                              table = TRUE)
cat("Weightr adjusted d:", weightr_result$adj.est, "\n")

# Step 6: Use puniform correction
puniform_result <- puniform(y = observed_d, v = se_d^2, method = "est")
cat("Punifom corrected d:", puniform_result$est, "\n")





## Scenario 1: Little bias

If the statistical power is high, there are few non-significant values, so the effect size distribution of only significant effects is not strongly inflated. Corrections like BUCSS and p-uniform lead to effect sizes close to the observed effect size. The grey line (statistical power) is high for effect sizes larger than X, and the corrected effect sizes remain close. When the observed power is 50%, the observed p-values are 0.05, and half of the findings will be non-significant in the long run.

```{r simulation 1}
#| echo: false
# 
# # Loop over mu1 from 0 to 1 in steps of 0.1
# mu1_values <- seq(0, 1, by = 0.01)
# 
# all_results <- lapply(mu1_values, function(m1) {
#   res <- compare_effect_size_corrections(
#     n = 300, 
#     mu1 = m1, 
#     mu2 = 0, 
#     sd = 1, 
#     sided = 2,
#     exact = TRUE,
#     retrodesign_d = NULL)
#   res$mu1 <- m1
#   res
# })
# 
# # Combine all into one table
# final_table <- bind_rows(all_results) %>%
#   select(mu1, everything())
# 
# # Reshape the data to long format
# plot_data <- final_table %>%
#   filter(significant) %>%  # Only include significant results
#   pivot_longer(
#     cols = c(power, BUCSS, puniform, type_m),
#     names_to = "Method",
#     values_to = "Adjusted_d"
#   )
# 
# # Plot
# ggplot(plot_data, aes(x = observed_d, y = Adjusted_d, color = Method)) +
#   #geom_point(aes(shape = Method), size = 2) +
#   geom_line(aes(group = Method, linetype = Method), linewidth = 2, show.legend = FALSE) +
#   scale_linetype_manual(values = c(
#   "solid", "solid", "dotted", "dotdash")) +
#   scale_color_manual(values = c("#E41A1C", "lightgrey", "black", "dodgerblue")) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   coord_fixed(ratio = 1) +
#   labs(
#     title = "Bias-Corrected Estimates vs. Observed d",
#     x = "Observed Cohen's d",
#     y = "Corrected Cohen's d and Power",
#     color = "Correction Method"
#   ) +
#   theme_minimal(base_size = 13)
# 

```

## Scenario 2: Lots of bias

If power is low overall, bias is high, and therefore the effect size adjustment d more considerable.

```{r simulation 2}
#| echo: false

# # Loop over mu1 from 0 to 1 in steps of 0.1
# mu1_values <- seq(0, 1, by = 0.01)
# 
# all_results <- lapply(mu1_values, function(m1) {
#   res <- compare_effect_size_corrections(
#     n = 40, 
#     mu1 = m1, 
#     mu2 = 0, 
#     sd = 1, 
#     sided = 2,
#     exact = TRUE,
#     retrodesign_d = NULL)
#   res$mu1 <- m1
#   res
# })
# 
# # Combine all into one table
# final_table <- bind_rows(all_results) %>%
#   select(mu1, everything())
# 
# # Recode method labels for consistency
# plot_data <- plot_data %>%
#   mutate(Method = recode(Method,
#                          "puniform" = "P-uniform",
#                          "type_m" = "Type M",
#                          "BUCSS" = "BUCSS",
#                          "power" = "Power"))
# 
# # Set factor levels to control legend and line order
# plot_data$Method <- factor(plot_data$Method,
#                            levels = c("Power", "Type M", "P-uniform", "BUCSS"))
# 
# # Plot
# ggplot(plot_data, aes(x = observed_d, y = Adjusted_d, color = Method, linetype = Method)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray80", linewidth = 1) +
#   geom_line(linewidth = 1.7) +
#   scale_color_manual(values = c(
#     "P-uniform" = "#bb7125",
#     "Type M" = "#ffdd00",
#     "Power" = "#00978d",
#     "BUCSS" = "#1c4286"
#   )) +
#   scale_linetype_manual(values = c(
#     "Power" = "solid",
#     "Type M" = "solid",
#     "P-uniform" = "solid",
#     "BUCSS" = "dotted"
#   )) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   coord_fixed(ratio = 1) +
#   labs(
#     title = "Bias-Corrected Estimates vs. Observed d",
#     x = "Cohen's d",
#     y = "Corrected Cohen's d and Power",
#     color = NULL,
#     linetype = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   labs(color = NULL, fill = NULL, linetype = NULL) + 
#   theme(
#     legend.position = "bottom",
#     text = element_text(family = "Segoe UI")
#   )
# 
# # Critical effect size
# # critical value from the t statistic
# m1<-0.5
# m2<-0
# sd1<-1
# sd2<-1
# n1<-40
# n2<-40
# 
# se <- sqrt(sd1^2 / n1 + sd2^2 / n2)
# t <- (m1 - m2) / se
# crit_res <- critical_t2s(t = t, n1 = n1, n2 = n2, se = se)
# crit_res$dc

```

## Scenario 3: Medium bias

The Type M error quantifies how much, on average, effect sizes are inflated. Gelman and Carlin [-@gelman2014] write in their example of a retrospective design analysis: "Conditional on the estimate being statistically significant, there is a 46% chance it will have the wrong sign (the Type S error rate), and in expectation the estimated effect will be 77 times too high (the exaggeration ratio). Both these conclusions can be reformulated more formally. First, the observed effect either has the wrong sign, or it hasn't, and in the long run effects selected for significance can be expected to have the wrong sign 46% of the time, and in the right direction 56% of the time. Furthermore, on average effect sizes will be 77 times to high, but the observed effect size is exaggerated by a different factor. After all, if the average absolute unbiased effect size estimate is 0.5, and the average bias effect size estimate is 1, the estimated effect is on average 2 times larger than the true effect size. In essence, the whole distribution needs to move 0.5 down to be unbiased. But an observed effect of 0.6 is moved down to 0.1, so this individual effect size is inflated 6 times from what it should be, and ab observed effect of 1.5 is moved down to 1, which is is an inflation ratio of 1.5. Where the average inflation can be used to educate researchers about the risk of inflation, it can not be used to adjust individual effect sizes for selection bias.

```{r simulation 3}
#| echo: false
# 
# # Loop over mu1 from 0 to 1 in steps of 0.1
# mu1_values <- seq(0, 1, by = 0.01)
# 
# all_results <- lapply(mu1_values, function(m1) {
#   res <- compare_effect_size_corrections(
#     n = 60, 
#     mu1 = m1, 
#     mu2 = 0, 
#     sd = 1, 
#     sided = 2,
#     exact = TRUE,
#     retrodesign_d = 0.5)
#   res$mu1 <- m1
#   res
# })
# 
# # Combine all into one table
# final_table <- bind_rows(all_results) %>%
#   select(mu1, everything())
# 
# # Reshape the data to long format
# plot_data <- final_table %>%
#   filter(significant) %>%  # Only include significant results
#   pivot_longer(
#     cols = c(power, BUCSS, puniform, type_m),
#     names_to = "Method",
#     values_to = "Adjusted_d"
#   )
# 
# # Plot
# ggplot(plot_data, aes(x = observed_d, y = Adjusted_d, color = Method)) +
#   #geom_point(aes(shape = Method), size = 2) +
#   geom_line(aes(group = Method, linetype = Method), linewidth = 2, show.legend = FALSE) +
#   scale_linetype_manual(values = c(
#   "solid", "solid", "dotted", "dotdash")) +
#   scale_color_manual(values = c("#E41A1C", "lightgrey", "black", "dodgerblue")) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   coord_fixed(ratio = 1) +
#   labs(
#     title = "Bias-Corrected Estimates vs. Observed d",
#     x = "Observed Cohen's d",
#     y = "Corrected Cohen's d and Power",
#     color = "Correction Method"
#   ) +
#   theme_minimal(base_size = 13)
# 

```

