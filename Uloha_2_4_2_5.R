# ==============================================================================
# ÚKOL 3 – ÚLOHY 2.4 a 2.5
# NKPC: NLS, ML a LR test dlouhodobé vertikality

# ------------------------------
# 1. Příprava dat pro 2.4 a 2.5
# ------------------------------
# Používáme stejnou datovou logiku jako v úlohách 2.2 a 2.3.
# Rozšířený hybridní model obsahuje pi_{t-2}, protože předchozí část řešila opomenutou dynamiku.

nkpc_data <- uk_data_model %>%
  mutate(
    pi      = inflation,
    pi_fwd  = lead(inflation, 1),
    pi_lag  = lag(inflation, 1),
    pi_lag2 = lag(inflation, 2),
    s_t     = ulc_gap
  ) %>%
  select(date, pi, pi_fwd, pi_lag, pi_lag2, s_t) %>%
  tidyr::drop_na()

cat("\nPočet pozorování pro NLS/ML:", nrow(nkpc_data), "\n")

# Matice pro rychlejší výpočty
Y  <- nkpc_data$pi
PF <- nkpc_data$pi_fwd
PL <- nkpc_data$pi_lag
P2 <- nkpc_data$pi_lag2
S  <- nkpc_data$s_t
n_obs <- length(Y)

# ------------------------------
# 2. Pomocné funkce
# ------------------------------
rmse <- function(resid) sqrt(mean(resid^2))

loglik_normal <- function(resid) {
  sigma <- sqrt(mean(resid^2))
  sum(dnorm(resid, mean = 0, sd = sigma, log = TRUE))
}

loglik_t <- function(resid, nu) {
  sigma <- sqrt(mean(resid^2))
  sum(dt(resid / sigma, df = nu, log = TRUE) - log(sigma))
}

lr_test <- function(ll_unrestr, ll_restr, df = 1) {
  LR <- 2 * (ll_unrestr - ll_restr)
  p  <- pchisq(LR, df = df, lower.tail = FALSE)
  c(LR = LR, p_value = p)
}

# Strukturální hybridní NKPC dle Galí-Gertler logiky:
# gamma_f = beta*theta / phi
# gamma_b = omega / phi
# lambda  = (1-omega)*(1-theta)*(1-beta*theta) / phi
# phi     = theta + omega*(1 - theta*(1-beta))
# Do rovnice přidáváme b_lag2*pi_{t-2} jako řešení opomenuté dynamiky z 2.2.
struct_coefs <- function(beta, theta, omega) {
  # optim() předává pojmenované vektory. Pokud se názvy neodstraní,
  # c(gamma_f = gamma_f) může vytvořit název typu "gamma_f.beta".
  # Pak cf["gamma_f"] vrací NA a optim hlásí: fn must be finite.
  beta  <- unname(as.numeric(beta))
  theta <- unname(as.numeric(theta))
  omega <- unname(as.numeric(omega))

  phi <- theta + omega * (1 - theta * (1 - beta))
  gamma_f <- beta * theta / phi
  gamma_b <- omega / phi
  lambda  <- (1 - omega) * (1 - theta) * (1 - beta * theta) / phi

  out <- c(gamma_f = gamma_f, gamma_b = gamma_b, lambda = lambda, phi = phi)
  if (any(!is.finite(out))) {
    out[] <- NA_real_
  }
  out
}

resid_reduced_unrestr <- function(par) {
  # par = bf, bb, kappa, b_lag2
  Y - par[1] * PF - par[2] * PL - par[3] * S - par[4] * P2
}

resid_reduced_restr <- function(par) {
  # Omezení dlouhodobé vertikality: bf + bb = 1 => bb = 1 - bf
  # par = bf, kappa, b_lag2
  bf <- par[1]
  bb <- 1 - bf
  Y - bf * PF - bb * PL - par[2] * S - par[3] * P2
}

resid_struct_unrestr <- function(par) {
  # par = beta, theta, omega, b_lag2
  beta  <- par[1]
  theta <- par[2]
  omega <- par[3]
  b2    <- par[4]
  cf <- struct_coefs(beta, theta, omega)
  Y - cf["gamma_f"] * PF - cf["gamma_b"] * PL - cf["lambda"] * S - b2 * P2
}

resid_struct_restr <- function(par) {
  # Omezení beta = 1; par = theta, omega, b_lag2
  beta  <- 1
  theta <- par[1]
  omega <- par[2]
  b2    <- par[3]
  cf <- struct_coefs(beta, theta, omega)
  Y - cf["gamma_f"] * PF - cf["gamma_b"] * PL - cf["lambda"] * S - b2 * P2
}

sse <- function(resid) {
  if (any(!is.finite(resid))) return(1e100)
  sum(resid^2)
}

# ------------------------------
# 3. Úloha 2.4 – NLS: redukovaná a strukturální podoba
# ------------------------------
cat("\n===============================================================================\n")
cat("2.4 NLS ODHADY\n")
cat("===============================================================================\n")

# 3.1 Redukovaná podoba – neomezená
nls_red_unrestr <- nls(
  pi ~ bf * pi_fwd + bb * pi_lag + kappa * s_t + b_lag2 * pi_lag2,
  data = nkpc_data,
  start = list(bf = 0.6, bb = 0.3, kappa = 0.1, b_lag2 = 0.1),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

# 3.2 Redukovaná podoba – omezená dlouhodobá vertikalita bf + bb = 1
nls_red_restr <- nls(
  pi ~ bf * pi_fwd + (1 - bf) * pi_lag + kappa * s_t + b_lag2 * pi_lag2,
  data = nkpc_data,
  start = list(bf = 0.6, kappa = 0.1, b_lag2 = 0.1),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

cat("\n--- NLS redukovaná podoba: neomezený model ---\n")
print(summary(nls_red_unrestr))
cat("\n--- NLS redukovaná podoba: omezený model bf + bb = 1 ---\n")
print(summary(nls_red_restr))

# 3.3 Strukturální podoba – NLS přes minimalizaci SSE s omezením parametrů
# beta povolujeme v intervalu (0.01, 1.50), theta a omega v (0.01, 0.99).
# To zabraňuje numerickým nesmyslům a zároveň umožňuje test beta = 1.

nls_struct_unrestr <- optim(
  par = c(beta = 0.99, theta = 0.70, omega = 0.30, b_lag2 = 0.05),
  fn = function(par) sse(resid_struct_unrestr(par)),
  method = "L-BFGS-B",
  lower = c(beta = 0.01, theta = 0.01, omega = 0.01, b_lag2 = -2),
  upper = c(beta = 1.50, theta = 0.99, omega = 0.99, b_lag2 =  2),
  control = list(maxit = 5000)
)

nls_struct_restr <- optim(
  par = c(theta = 0.70, omega = 0.30, b_lag2 = 0.05),
  fn = function(par) sse(resid_struct_restr(par)),
  method = "L-BFGS-B",
  lower = c(theta = 0.01, omega = 0.01, b_lag2 = -2),
  upper = c(theta = 0.99, omega = 0.99, b_lag2 =  2),
  control = list(maxit = 5000)
)

cat("\n--- NLS strukturální podoba: neomezený model ---\n")
print(nls_struct_unrestr$par)
print(struct_coefs(nls_struct_unrestr$par["beta"],
                   nls_struct_unrestr$par["theta"],
                   nls_struct_unrestr$par["omega"]))
cat("SSE:", nls_struct_unrestr$value, " RMSE:", rmse(resid_struct_unrestr(nls_struct_unrestr$par)), "\n")

cat("\n--- NLS strukturální podoba: omezený model beta = 1 ---\n")
print(c(beta = 1, nls_struct_restr$par))
print(struct_coefs(1,
                   nls_struct_restr$par["theta"],
                   nls_struct_restr$par["omega"]))
cat("SSE:", nls_struct_restr$value, " RMSE:", rmse(resid_struct_restr(nls_struct_restr$par)), "\n")

# ------------------------------
# 4. Úloha 2.4 – ML: normální a t-rozdělené chyby
# ------------------------------
cat("\n===============================================================================\n")
cat("2.4 ML ODHADY – NORMÁLNÍ A t-ROZDĚLENÉ CHYBY\n")
cat("===============================================================================\n")

# 4.1 ML redukovaná podoba – normální chyby
ml_red_norm_unrestr <- optim(
  par = c(coef(nls_red_unrestr), log_sigma = log(sd(resid(nls_red_unrestr)))),
  fn = function(par) {
    r <- resid_reduced_unrestr(par[1:4])
    sigma <- exp(par[5])
    val <- -sum(dnorm(r, mean = 0, sd = sigma, log = TRUE))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(-5, -5, -5, -5, -10),
  upper = c( 5,  5,  5,  5,  10),
  control = list(maxit = 5000)
)

ml_red_norm_restr <- optim(
  par = c(coef(nls_red_restr), log_sigma = log(sd(resid(nls_red_restr)))),
  fn = function(par) {
    r <- resid_reduced_restr(par[1:3])
    sigma <- exp(par[4])
    val <- -sum(dnorm(r, mean = 0, sd = sigma, log = TRUE))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(-5, -5, -5, -10),
  upper = c( 5,  5,  5,  10),
  control = list(maxit = 5000)
)

# 4.2 ML redukovaná podoba – t chyby
ml_red_t_unrestr <- optim(
  par = c(coef(nls_red_unrestr), log_sigma = log(sd(resid(nls_red_unrestr))), log_nu_minus2 = log(8)),
  fn = function(par) {
    r <- resid_reduced_unrestr(par[1:4])
    sigma <- exp(par[5])
    nu <- 2 + exp(par[6])
    val <- -sum(dt(r / sigma, df = nu, log = TRUE) - log(sigma))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(-5, -5, -5, -5, -10, log(0.1)),
  upper = c( 5,  5,  5,  5,  10, log(200)),
  control = list(maxit = 5000)
)

ml_red_t_restr <- optim(
  par = c(coef(nls_red_restr), log_sigma = log(sd(resid(nls_red_restr))), log_nu_minus2 = log(8)),
  fn = function(par) {
    r <- resid_reduced_restr(par[1:3])
    sigma <- exp(par[4])
    nu <- 2 + exp(par[5])
    val <- -sum(dt(r / sigma, df = nu, log = TRUE) - log(sigma))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(-5, -5, -5, -10, log(0.1)),
  upper = c( 5,  5,  5,  10, log(200)),
  control = list(maxit = 5000)
)

# 4.3 ML strukturální podoba – normální chyby
ml_struct_norm_unrestr <- optim(
  par = c(nls_struct_unrestr$par, log_sigma = log(sd(resid_struct_unrestr(nls_struct_unrestr$par)))),
  fn = function(par) {
    r <- resid_struct_unrestr(par[1:4])
    sigma <- exp(par[5])
    val <- -sum(dnorm(r, mean = 0, sd = sigma, log = TRUE))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, 0.01, -2, -10),
  upper = c(1.50, 0.99, 0.99,  2,  10),
  control = list(maxit = 5000)
)

ml_struct_norm_restr <- optim(
  par = c(nls_struct_restr$par, log_sigma = log(sd(resid_struct_restr(nls_struct_restr$par)))),
  fn = function(par) {
    r <- resid_struct_restr(par[1:3])
    sigma <- exp(par[4])
    val <- -sum(dnorm(r, mean = 0, sd = sigma, log = TRUE))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, -2, -10),
  upper = c(0.99, 0.99,  2,  10),
  control = list(maxit = 5000)
)

# 4.4 ML strukturální podoba – t chyby
ml_struct_t_unrestr <- optim(
  par = c(nls_struct_unrestr$par,
          log_sigma = log(sd(resid_struct_unrestr(nls_struct_unrestr$par))),
          log_nu_minus2 = log(8)),
  fn = function(par) {
    r <- resid_struct_unrestr(par[1:4])
    sigma <- exp(par[5])
    nu <- 2 + exp(par[6])
    val <- -sum(dt(r / sigma, df = nu, log = TRUE) - log(sigma))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, 0.01, -2, -10, log(0.1)),
  upper = c(1.50, 0.99, 0.99,  2,  10, log(200)),
  control = list(maxit = 5000)
)

ml_struct_t_restr <- optim(
  par = c(nls_struct_restr$par,
          log_sigma = log(sd(resid_struct_restr(nls_struct_restr$par))),
          log_nu_minus2 = log(8)),
  fn = function(par) {
    r <- resid_struct_restr(par[1:3])
    sigma <- exp(par[4])
    nu <- 2 + exp(par[5])
    val <- -sum(dt(r / sigma, df = nu, log = TRUE) - log(sigma))
    if (!is.finite(val)) 1e100 else val
  },
  method = "L-BFGS-B",
  lower = c(0.01, 0.01, -2, -10, log(0.1)),
  upper = c(0.99, 0.99,  2,  10, log(200)),
  control = list(maxit = 5000)
)

# Výpis ML parametrů
cat("\n--- ML redukovaná podoba, normální chyby, neomezený model ---\n")
print(ml_red_norm_unrestr$par)
cat("LogLik:", -ml_red_norm_unrestr$value, "\n")

cat("\n--- ML redukovaná podoba, t chyby, neomezený model ---\n")
print(c(ml_red_t_unrestr$par[1:5], nu = 2 + exp(ml_red_t_unrestr$par[6])))
cat("LogLik:", -ml_red_t_unrestr$value, "\n")

cat("\n--- ML strukturální podoba, normální chyby, neomezený model ---\n")
print(ml_struct_norm_unrestr$par)
cat("LogLik:", -ml_struct_norm_unrestr$value, "\n")

cat("\n--- ML strukturální podoba, t chyby, neomezený model ---\n")
print(c(ml_struct_t_unrestr$par[1:5], nu = 2 + exp(ml_struct_t_unrestr$par[6])))
cat("LogLik:", -ml_struct_t_unrestr$value, "\n")


round_df <- function(df, digits = 4) {
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  df
}

sig_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", ""))))
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits, decimal.mark = ","))
}

fmt_est <- function(est, p = NA_real_, digits = 4) {
  paste0(fmt_num(est, digits), sig_stars(p))
}

extract_gmm_pvalues <- function(gmm_obj, par_names) {
  out <- setNames(rep(NA_real_, length(par_names)), par_names)
  sm <- tryCatch(summary(gmm_obj), error = function(e) NULL)
  if (is.null(sm)) return(out)

  cm <- tryCatch(coef(sm), error = function(e) NULL)
  if (is.null(cm)) cm <- tryCatch(sm$coefficients, error = function(e) NULL)
  if (is.null(cm) || is.null(dim(cm))) return(out)

  p_col <- grep("Pr|p.value|p-value|P", colnames(cm), ignore.case = TRUE, value = TRUE)
  if (length(p_col) == 0) return(out)
  p_col <- p_col[1]

  common <- intersect(par_names, rownames(cm))
  out[common] <- cm[common, p_col]
  out
}

# ------------------------------
# 5. Úloha 2.5 – srovnání metod a LR testy
# ------------------------------
cat("\n===============================================================================\n")
cat("2.5 SROVNÁNÍ NLS, ML A GMM + LR TESTY DLOUHODOBÉ VERTIKALITY\n")
cat("===============================================================================\n")

# 5.1 Tabulka redukované podoby
red_table <- data.frame(
  Metoda = c("GMM rozšířený hybrid", "NLS", "ML normální", "ML t"),
  bf = NA_real_, bb = NA_real_, kappa = NA_real_, b_lag2 = NA_real_,
  beta_longrun = NA_real_, RMSE = NA_real_, LogLik = NA_real_, nu_t = NA_real_
)

if (exists("gmm_hyb_ext") && exists("gmm_data_mat2")) {
  cg <- coef(gmm_hyb_ext)
  rg <- gmm_data_mat2[, "pi"] -
    cg["b_f"]    * gmm_data_mat2[, "pi_fwd"] -
    cg["b_b"]    * gmm_data_mat2[, "pi_lag"] -
    cg["kappa"]  * gmm_data_mat2[, "s_t"] -
    cg["b_lag2"] * gmm_data_mat2[, "pi_lag2"]
  red_table[1, c("bf", "bb", "kappa", "b_lag2")] <- c(cg["b_f"], cg["b_b"], cg["kappa"], cg["b_lag2"])
  red_table[1, "beta_longrun"] <- cg["b_f"] + cg["b_b"]
  red_table[1, "RMSE"] <- rmse(rg)
}

cn <- coef(nls_red_unrestr)
rn <- resid_reduced_unrestr(c(cn["bf"], cn["bb"], cn["kappa"], cn["b_lag2"]))
red_table[2, c("bf", "bb", "kappa", "b_lag2")] <- c(cn["bf"], cn["bb"], cn["kappa"], cn["b_lag2"])
red_table[2, "beta_longrun"] <- cn["bf"] + cn["bb"]
red_table[2, "RMSE"] <- rmse(rn)
red_table[2, "LogLik"] <- loglik_normal(rn)

pm <- ml_red_norm_unrestr$par
rm <- resid_reduced_unrestr(pm[1:4])
red_table[3, c("bf", "bb", "kappa", "b_lag2")] <- pm[1:4]
red_table[3, "beta_longrun"] <- pm[1] + pm[2]
red_table[3, "RMSE"] <- rmse(rm)
red_table[3, "LogLik"] <- -ml_red_norm_unrestr$value

pt <- ml_red_t_unrestr$par
rt <- resid_reduced_unrestr(pt[1:4])
red_table[4, c("bf", "bb", "kappa", "b_lag2")] <- pt[1:4]
red_table[4, "beta_longrun"] <- pt[1] + pt[2]
red_table[4, "RMSE"] <- rmse(rt)
red_table[4, "LogLik"] <- -ml_red_t_unrestr$value
red_table[4, "nu_t"] <- 2 + exp(pt[6])

cat("\n--- Srovnání redukované podoby ---\n")
print(round_df(red_table, 4))

# 5.1b Tabulka NLS redukované podoby s hvězdičkami významnosti
# Tuto tabulku použijte v textu tam, kde chcete ukázat statistickou významnost NLS odhadu.
nls_red_unrestr_sum <- coef(summary(nls_red_unrestr))
nls_red_restr_sum   <- coef(summary(nls_red_restr))

nls_red_table_stars <- data.frame(
  Parametr = c("beta_f", "beta_b", "kappa", "beta_lag2", "beta_f + beta_b", "RMSE"),
  Neomezeny_model = c(
    fmt_est(nls_red_unrestr_sum["bf", "Estimate"], nls_red_unrestr_sum["bf", "Pr(>|t|)"]),
    fmt_est(nls_red_unrestr_sum["bb", "Estimate"], nls_red_unrestr_sum["bb", "Pr(>|t|)"]),
    fmt_est(nls_red_unrestr_sum["kappa", "Estimate"], nls_red_unrestr_sum["kappa", "Pr(>|t|)"]),
    fmt_est(nls_red_unrestr_sum["b_lag2", "Estimate"], nls_red_unrestr_sum["b_lag2", "Pr(>|t|)"]),
    fmt_num(sum(coef(nls_red_unrestr)[c("bf", "bb")]), 4),
    fmt_num(rmse(resid(nls_red_unrestr)), 4)
  ),
  Omezeny_model = c(
    fmt_est(nls_red_restr_sum["bf", "Estimate"], nls_red_restr_sum["bf", "Pr(>|t|)"]),
    paste0(fmt_num(1 - coef(nls_red_restr)["bf"], 4), " (impl.)"),
    fmt_est(nls_red_restr_sum["kappa", "Estimate"], nls_red_restr_sum["kappa", "Pr(>|t|)"]),
    fmt_est(nls_red_restr_sum["b_lag2", "Estimate"], nls_red_restr_sum["b_lag2", "Pr(>|t|)"]),
    "1,0000",
    fmt_num(rmse(resid(nls_red_restr)), 4)
  ),
  stringsAsFactors = FALSE
)

cat("\n--- NLS redukovaná podoba s hvězdičkami významnosti ---\n")
print(nls_red_table_stars, row.names = FALSE)
cat("Pozn.: V omezeném modelu je beta_b dopočteno jako 1 - beta_f, proto není samostatně testováno.\n")
cat("Pozn.: *** p<0,01; ** p<0,05; * p<0,10.\n")

# 5.1c Srovnávací tabulka s hvězdičkami tam, kde jsou dostupné p-hodnoty
# GMM: hvězdičky pouze pokud summary(gmm_hyb_ext) poskytne p-hodnoty.
# NLS: hvězdičky ze summary(nls). ML řádky necháváme bez hvězdiček, protože zde nepočítáme SE.
red_table_stars <- red_table
red_table_stars[] <- lapply(red_table_stars, as.character)

if (exists("gmm_hyb_ext")) {
  gmm_p <- extract_gmm_pvalues(gmm_hyb_ext, c("b_f", "b_b", "kappa", "b_lag2"))
  red_table_stars[1, "bf"]     <- fmt_est(red_table[1, "bf"],     gmm_p["b_f"])
  red_table_stars[1, "bb"]     <- fmt_est(red_table[1, "bb"],     gmm_p["b_b"])
  red_table_stars[1, "kappa"]  <- fmt_est(red_table[1, "kappa"],  gmm_p["kappa"])
  red_table_stars[1, "b_lag2"] <- fmt_est(red_table[1, "b_lag2"], gmm_p["b_lag2"])
}

red_table_stars[2, "bf"]     <- fmt_est(red_table[2, "bf"],     nls_red_unrestr_sum["bf", "Pr(>|t|)"])
red_table_stars[2, "bb"]     <- fmt_est(red_table[2, "bb"],     nls_red_unrestr_sum["bb", "Pr(>|t|)"])
red_table_stars[2, "kappa"]  <- fmt_est(red_table[2, "kappa"],  nls_red_unrestr_sum["kappa", "Pr(>|t|)"])
red_table_stars[2, "b_lag2"] <- fmt_est(red_table[2, "b_lag2"], nls_red_unrestr_sum["b_lag2", "Pr(>|t|)"])

# Formátování ostatních numerických sloupců bez hvězdiček
for (j in c("beta_longrun", "RMSE", "LogLik", "nu_t")) {
  red_table_stars[[j]] <- ifelse(is.na(red_table[[j]]), "--", fmt_num(red_table[[j]], 4))
}
for (i in c(3, 4)) {
  for (j in c("bf", "bb", "kappa", "b_lag2")) {
    red_table_stars[i, j] <- fmt_num(red_table[i, j], 4)
  }
}

cat("\n--- Srovnání redukované podoby s hvězdičkami tam, kde jsou dostupné p-hodnoty ---\n")
print(red_table_stars, row.names = FALSE)
cat("Pozn.: Hvězdičky jsou uvedeny jen pro GMM/NLS, pokud jsou dostupné p-hodnoty. ML řádky jsou bez hvězdiček.\n")
cat("Pozn.: *** p<0,01; ** p<0,05; * p<0,10.\n")


# 5.2 Tabulka strukturální podoby
struct_table <- data.frame(
  Metoda = c("NLS", "ML normální", "ML t"),
  beta = NA_real_, theta = NA_real_, omega = NA_real_,
  gamma_f = NA_real_, gamma_b = NA_real_, lambda = NA_real_,
  b_lag2 = NA_real_, beta_longrun_reduced = NA_real_, RMSE = NA_real_, LogLik = NA_real_, nu_t = NA_real_
)

ps <- nls_struct_unrestr$par
cf <- struct_coefs(ps["beta"], ps["theta"], ps["omega"])
rs <- resid_struct_unrestr(ps)
struct_table[1, c("beta", "theta", "omega", "b_lag2")] <- c(ps["beta"], ps["theta"], ps["omega"], ps["b_lag2"])
struct_table[1, c("gamma_f", "gamma_b", "lambda")] <- cf[c("gamma_f", "gamma_b", "lambda")]
struct_table[1, "beta_longrun_reduced"] <- cf["gamma_f"] + cf["gamma_b"]
struct_table[1, "RMSE"] <- rmse(rs)
struct_table[1, "LogLik"] <- loglik_normal(rs)

psm <- ml_struct_norm_unrestr$par
cfm <- struct_coefs(psm[1], psm[2], psm[3])
rsm <- resid_struct_unrestr(psm[1:4])
struct_table[2, c("beta", "theta", "omega", "b_lag2")] <- psm[1:4]
struct_table[2, c("gamma_f", "gamma_b", "lambda")] <- cfm[c("gamma_f", "gamma_b", "lambda")]
struct_table[2, "beta_longrun_reduced"] <- cfm["gamma_f"] + cfm["gamma_b"]
struct_table[2, "RMSE"] <- rmse(rsm)
struct_table[2, "LogLik"] <- -ml_struct_norm_unrestr$value

pst <- ml_struct_t_unrestr$par
cft <- struct_coefs(pst[1], pst[2], pst[3])
rst <- resid_struct_unrestr(pst[1:4])
struct_table[3, c("beta", "theta", "omega", "b_lag2")] <- pst[1:4]
struct_table[3, c("gamma_f", "gamma_b", "lambda")] <- cft[c("gamma_f", "gamma_b", "lambda")]
struct_table[3, "beta_longrun_reduced"] <- cft["gamma_f"] + cft["gamma_b"]
struct_table[3, "RMSE"] <- rmse(rst)
struct_table[3, "LogLik"] <- -ml_struct_t_unrestr$value
struct_table[3, "nu_t"] <- 2 + exp(pst[6])

cat("\n--- Srovnání strukturální podoby ---\n")
print(round_df(struct_table, 4))

# 5.3 LR testy dlouhodobé vertikality
# Redukovaná podoba: H0 bf + bb = 1.
# Strukturální podoba: H0 beta = 1.

ll_nls_red_unrestr <- loglik_normal(resid_reduced_unrestr(coef(nls_red_unrestr)))
ll_nls_red_restr   <- loglik_normal(resid_reduced_restr(coef(nls_red_restr)))
ll_nls_struct_unrestr <- loglik_normal(resid_struct_unrestr(nls_struct_unrestr$par))
ll_nls_struct_restr   <- loglik_normal(resid_struct_restr(nls_struct_restr$par))

lr_results <- data.frame(
  Model = c(
    "NLS redukovaná: H0 bf+bb=1",
    "NLS strukturální: H0 beta=1",
    "ML normální redukovaná: H0 bf+bb=1",
    "ML normální strukturální: H0 beta=1",
    "ML t redukovaná: H0 bf+bb=1",
    "ML t strukturální: H0 beta=1"
  ),
  LR = NA_real_,
  p_value = NA_real_
)

lr_results[1, c("LR", "p_value")] <- lr_test(ll_nls_red_unrestr, ll_nls_red_restr, df = 1)
lr_results[2, c("LR", "p_value")] <- lr_test(ll_nls_struct_unrestr, ll_nls_struct_restr, df = 1)
lr_results[3, c("LR", "p_value")] <- lr_test(-ml_red_norm_unrestr$value, -ml_red_norm_restr$value, df = 1)
lr_results[4, c("LR", "p_value")] <- lr_test(-ml_struct_norm_unrestr$value, -ml_struct_norm_restr$value, df = 1)
lr_results[5, c("LR", "p_value")] <- lr_test(-ml_red_t_unrestr$value, -ml_red_t_restr$value, df = 1)
lr_results[6, c("LR", "p_value")] <- lr_test(-ml_struct_t_unrestr$value, -ml_struct_t_restr$value, df = 1)

cat("\n--- LR testy dlouhodobé vertikality ---\n")
print(round_df(lr_results, 4))

