
gum_model1 <- list(
    type = "gumbel",
    npar = 2,
    alpha_ind = "alpha_ind <- 2",
    lambda1 = "lambda <-  beta[1] * df$points",
    lambda2 = "lambda <-  beta[1] * df$points"
)

gum_model2 <- list(
    type = "gumbel",
    npar = 3,
    alpha_ind = "alpha_ind <- 3",
    lambda1 = "lambda <-   beta[1] * df$points  + beta[2] * df$bib1/nrow(df)",
    lambda2 = "lambda <-   beta[1] * df$points  + beta[2] * df$bib2/nrow(df)"
)

gum_model3 <- list(
    type = "gumbel",
    npar = "n_races + 2",
    alpha_ind = "alpha_ind <- n_races + 2",
    lambda1 = "lambda <-   beta[1] * df$points  + beta[race_number + 1] * df$bib1/nrow(df)",
    lambda2 = "lambda <-   beta[1] * df$points  + beta[race_number + 1] * df$bib2/nrow(df)"
)

gum_model4 <- list(
    type = "gumbel",
    npar = "n_races + 3",
    alpha_ind = "alpha_ind <- n_races + 3",
    lambda1 = "lambda <-   beta[1] * df$points  + beta[race_number + 2] * df$bib1/nrow(df)",
    lambda2 = "lambda <-   beta[2] * df$points  + beta[race_number + 2] * df$bib2/nrow(df)"
)

gum_model5 <- list(
    type = "gumbel",
    npar = "n_races + 1",
    alpha_ind = "alpha_ind <- n_races + 1",
    lambda1 = "lambda <-   beta[race_number] * df$bib1/nrow(df)",
    lambda2 = "lambda <-   beta[race_number] * df$bib2/nrow(df)"
)

gum_model6 <- list(
    type = "gumbel",
    npar = "n_races + 2",
    alpha_ind = "alpha_ind <- n_races + 2",
    lambda1 = "lambda <-   beta[1] * df$points  + beta[race_number + 1] * df$bib1/nrow(df)",
    lambda2 = "lambda <-   beta[race_number + 1] * df$bib2/nrow(df)"
)

exp_model1 <- list(
    type = "exp",
    npar = 2,
    lambda1 = "lambda <-  exp(beta[1] + beta[2] * df$points)",
    lambda2 = "lambda <-  exp(beta[1] + beta[2] * df$points)"
)

exp_model2 <- list(
    type = "exp",
    npar = 3,
    lambda1 = "lambda <-  exp(beta[1] + beta[2] * df$points + beta[3] * df$bib1/nrow(df))",
    lambda2 = "lambda <-  exp(beta[1] + beta[2] * df$points + beta[3] * df$bib2/nrow(df))"
)

exp_model3 <- list(
    type = "exp",
    npar = "n_races + 2",
    lambda1 = "lambda <-   exp(beta[1] + beta[2] * df$points  + beta[race_number + 2] * df$bib1/nrow(df))",
    lambda2 = "lambda <-   exp(beta[1] + beta[2] * df$points  + beta[race_number + 2] * df$bib2/nrow(df))"
)

exp_model4 <- list(
    type = "exp",
    npar = "n_races + 4",
    lambda1 = "lambda <-   exp(beta[1] + beta[3] * df$points  + beta[race_number + 4] * df$bib1/nrow(df))",
    lambda2 = "lambda <-   exp(beta[2] + beta[4] * df$points  + beta[race_number + 4] * df$bib2/nrow(df))"
)

exp_model5 <- list(
    type = "exp",
    npar = "n_races + 2",
    lambda1 = "lambda <-   exp(beta[1] + beta[race_number + 2] * df$bib1/nrow(df))",
    lambda2 = "lambda <-   exp(beta[2] + beta[race_number + 2] * df$bib2/nrow(df))"
)

exp_model6 <- list(
    type = "exp",
    npar = "n_races + 3",
    lambda1 = "lambda <-   exp(beta[1] + beta[3] * df$points  + beta[race_number + 3] * df$bib1/nrow(df))",
    lambda2 = "lambda <-   exp(beta[2] + beta[race_number + 3] * df$bib2/nrow(df))"
)
