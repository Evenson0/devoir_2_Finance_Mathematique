# ********************************** Devoir 2 **********************************

# **************************** Finance mathématique ****************************

# Auteurs : Jad Ramy - Evenson Auguste
# Date : 12 Décembre 2024

#Exercice 2
# a)

BSOptionPrice <- function(S, K, r, T_t, sigma, isput) {
  # Calcul des paramètres d1 et d2
  d1 <- (log(S / K) + (r + (sigma^2) / 2) * T_t) / (sigma * sqrt(T_t))
  d2 <- d1 - sigma * sqrt(T_t)
  
  # Calcul du prix de l'option
  if (isput) {
    # Option de vente (put)
    option_price <- K * exp(-r * T_t) * pnorm(-d2) - S * pnorm(-d1)
  } else {
    # Option d'achat (call)
    option_price <- S * pnorm(d1) - K * exp(-r * T_t) * pnorm(d2)
  }
  
  return(option_price)
}

# b)

#option d'achat
option_price <- BSOptionPrice(S = 100, K = 105, r = 0.02, T_t = 0.5, sigma = 0.2, isput = FALSE)
cat("Le prix de l'option d'achat est :", option_price, "\n")
#option de vente
option_price <- BSOptionPrice(S = 100, K = 105, r = 0.02, T_t = 0.5, sigma = 0.2, isput = TRUE)
cat("Le prix de l'option de vente est :", option_price, "\n")

# c)

# Fonction pour trouver la volatilité implicite
BSImplicitVol <- function(OptionPrice, S, K, r, T_t, isput, tol = 1e-6) {
  # Fonction à minimiser : différence entre le prix observé et le prix Black-Scholes
  option_price_diff <- function(sigma) {
    BSOptionPrice(S, K, r, T_t, sigma, isput) - OptionPrice
  }
  
  # Utiliser uniroot pour trouver la volatilité implicite
  result <- tryCatch(
    {
      uniroot(option_price_diff, lower = 1e-6, upper = 5, tol = tol)
    },
    error = function(e) {
      return(NULL)  # Retourne NULL si aucune solution n'est trouvée
    }
  )
  
  if (!is.null(result)) {
    implied_vol <- result$root
    print(paste("La volatilité implicite est :", round(implied_vol, 6)))  # Affiche directement la volatilité implicite
    return(implied_vol)  # Retourne la volatilité implicite
  } else {
    print("Impossible de trouver une volatilité implicite.")
    return(NA)  # Retourne NA si uniroot échoue
  }
}

# d)

call_vol <- BSImplicitVol(OptionPrice = 2.7852, S = 100, K = 105, r = 0.02, T_t = 0.25, isput = FALSE)

# e)

put_vol <- BSImplicitVol(OptionPrice = 6.8249, S = 100, K = 105, r = 0.02, T_t = 0.75, isput = TRUE)

# Exercice 3

# a)

BinOptionPrice <- function(S, K, r, T_t, mu, sigma, n, isput) {
  # Calcul des paramètres
  h <- T_t / n
  u <- exp(mu * h + sigma * sqrt(h))
  d <- exp(mu * h - sigma * sqrt(h))
  p <- (exp(r * h) - d) / (u - d)
  disc <- exp(-r * h)
  
  # Prix de l'actif à maturité (t = T)
  ST <- numeric(n+1)
  for (j in 0:n) {
    ST[j+1] <- S * (u^j) * (d^(n-j))
  }
  
  # Payoffs à maturité
  payoff <- if (!isput) {
    # Call
    pmax(ST - K, 0)
  } else {
    # Put
    pmax(K - ST, 0)
  }
  
  # Roll-back vers t=0
  for (i in (n-1):0) {
    for (j in 0:i) {
      payoff[j+1] <- disc * (p * payoff[j+2] + (1 - p) * payoff[j+1])
    }
  }
  
  return(payoff[1])
}

# b)

S0 <- 100
K <- 105
r <- 0.02
T <- 0.5
sigma <- 0.2
mu <- r - sigma^2/2
n <- 20

# Prix d'un call
call_price <- BinOptionPrice(S0, K, r, T, mu, sigma, n, isput = FALSE)

# Prix d'un put
put_price <- BinOptionPrice(S0, K, r, T, mu, sigma, n, isput = TRUE)

call_price
put_price

# c)

# Paramètres fixes
S0 <- 100
K <- 105
r <- 0.02
T_t <- 0.5
sigma <- 0.2
isput <- FALSE  # Option d'achat (call)

# Valeurs spécifiques pour mu et n
mu_values <- c(-0.025, -0.01, 0, 0.01, 0.02, 0.025, 0.03, 0.05, 0.25, 0.5)
n_values <- c(10, 50, 100, 500, 1000, 2000, 5000, 10000)

# Prix Black-Scholes
bs_price <- BSOptionPrice(S0, K, r, T_t, sigma, isput)

# Initialisation des résultats
results <- data.frame(
  mu = numeric(),
  n = numeric(),
  `Binomial Price` = numeric(),
  `Black-Scholes Price` = numeric(),
  `Difference` = numeric()
)

# Boucles pour les combinaisons de mu et n
for (mu in mu_values) {
  for (n in n_values) {
    bin_price <- BinOptionPrice(S0, K, r, T_t, mu, sigma, n, isput)
    diff_from_bs <- abs(bin_price - bs_price)
    results <- rbind(
      results,
      data.frame(
        mu = mu,
        n = n,
        `Binomial Price` = bin_price,
        `Black-Scholes Price` = bs_price,
        `Difference` = diff_from_bs
      )
    )
  }
}

# Afficher les résultats
print(results)

# e)

# Paramètres
S0 <- 100
K <- 105
r <- 0.02
T <- 0.5
sigma <- 0.2
n <- 1000  # Nombre de pas dans l'arbre binomial

# Construction de l'arbre binomial (neutre au risque)
h <- T / n
u <- exp(sigma * sqrt(h))
d <- 1 / u
p <- (exp(r * h) - d) / (u - d)

# Calcul des valeurs possibles de ST et leurs probabilités
ST_values <- sapply(0:n, function(j) S0 * (u^j) * (d^(n - j)))
ST_probs <- sapply(0:n, function(j) choose(n, j) * (p^j) * ((1 - p)^(n - j)))

# Distribution lognormale
mu <- log(S0) + (r - sigma^2 / 2) * T
variance <- sigma^2 * T
lognorm_density <- dlnorm(ST_values, meanlog = mu, sdlog = sqrt(variance))

# Création du tableau
results <- data.frame(
  ST = ST_values,
  `Binomial Density` = ST_probs,
  `Lognormal Density` = lognorm_density
)

# Affichage du tableau
print(results)

# Exercice 4

# a)

SimulateStockPaths <- function(S, mu, r, sigma, nsteps, sizesteps, npaths, theseed, isPmeasure) {
  # Fixer la graine aléatoire pour la reproductibilité
  set.seed(theseed)
  
  # Initialisation
  dt <- sizesteps                     # Taille d'un pas de temps
  N <- nsteps                         # Nombre de pas de temps
  M <- npaths                         # Nombre de trajectoires
  paths <- matrix(0, nrow = M, ncol = N + 1)  # Matrice pour stocker les trajectoires
  paths[, 1] <- S                     # Initialisation avec le prix initial S0
  
  # Définir le drift en fonction de la mesure
  drift <- if (isPmeasure) {
    (mu - 0.5 * sigma^2) * dt
  } else {
    (r - 0.5 * sigma^2) * dt
  }
  
  # Simulation des trajectoires
  for (i in 1:N) {
    Z <- rnorm(M)                     # Simulation des termes aléatoires Z ~ N(0, 1)
    diffusion <- sigma * sqrt(dt) * Z
    paths[, i + 1] <- paths[, i] * exp(drift + diffusion)
  }
  
  return(paths)
}

# b)

SimulateStockPathsAntithetic <- function(S, mu, r, sigma, nsteps, sizesteps, npaths, theseed, isPmeasure) {
  # Vérification que npaths est pair
  if (npaths %% 2 != 0) {
    stop("Le nombre de trajectoires (npaths) doit être pair pour utiliser la méthode antithétique.")
  }
  
  # Fixer la graine pour la reproductibilité
  set.seed(theseed)
  
  # Initialisation
  dt <- sizesteps
  N <- nsteps
  M <- npaths
  halfM <- M / 2
  
  # Matrice des trajectoires
  paths <- matrix(0, nrow = M, ncol = N + 1)
  paths[, 1] <- S
  
  # Définir le drift
  drift <- if (isPmeasure) {
    (mu - 0.5 * sigma^2) * dt
  } else {
    (r - 0.5 * sigma^2) * dt
  }
  
  # Simulation des trajectoires standard pour M/2 premières trajectoires
  for (i in 1:N) {
    Z <- rnorm(halfM)
    diffusion <- sigma * sqrt(dt) * Z
    paths[1:halfM, i + 1] <- paths[1:halfM, i] * exp(drift + diffusion)
  }
  
  
  # Ré-initialisons le tout, pour clarté :
  set.seed(theseed)
  paths <- matrix(0, nrow = M, ncol = N + 1)
  paths[, 1] <- S
  
  Z_matrix <- matrix(0, nrow = halfM, ncol = N) # Pour stocker les Z
  # Première boucle : trajectoires normales
  for (i in 1:N) {
    Z <- rnorm(halfM)
    Z_matrix[, i] <- Z
    diffusion <- sigma * sqrt(dt) * Z
    paths[1:halfM, i + 1] <- paths[1:halfM, i] * exp(drift + diffusion)
  }
  
  # Deuxième boucle : trajectoires antithétiques (en utilisant -Z)
  for (i in 1:N) {
    Z_antithetic <- -Z_matrix[, i]
    diffusion_antithetic <- sigma * sqrt(dt) * Z_antithetic
    paths[(halfM+1):M, i + 1] <- paths[(halfM+1):M, i] * exp(drift + diffusion_antithetic)
  }
  
  return(paths)
}

# Exemple d'utilisation
S <- 100          # Prix initial
mu <- 0.05        # Drift sous la mesure réelle
r <- 0.03         # Taux sans risque
sigma <- 0.2      # Volatilité
nsteps <- 100     # Nombre de pas de temps
sizesteps <- 0.01 # Taille d'un pas de temps
npaths <- 1000    # Nombre de trajectoires (doit être pair)
theseed <- 42      # Graine
isPmeasure <- TRUE # Simulation sous la mesure P

# Simulation avec antithétiques
trajectories_antithetic <- SimulateStockPathsAntithetic(S, mu, r, sigma, nsteps, sizesteps, npaths, theseed, isPmeasure)

# Affichage des premières lignes
head(trajectories_antithetic)

#c)

# Paramètres
S0 <- 100
mu <- 0.07
r <- 0.02
sigma <- 0.20
T <- 0.5
dt <- 1/52
nsteps <- T/dt  # 0.5 * 52 = 26
npaths <- 10000
theseed <- 12345 # Graine choisie

# Simulation sous la mesure P (réelle)
StockPathsP <- SimulateStockPaths(
  S = S0, mu = mu, r = r, sigma = sigma,
  nsteps = nsteps, sizesteps = dt, npaths = npaths,
  theseed = theseed, isPmeasure = TRUE
)

StockPathsPantithetic <- SimulateStockPathsAntithetic(
  S = S0, mu = mu, r = r, sigma = sigma,
  nsteps = nsteps, sizesteps = dt, npaths = npaths,
  theseed = theseed, isPmeasure = TRUE
)

# Simulation sous la mesure Q (risque-neutre)
StockPathsQ <- SimulateStockPaths(
  S = S0, mu = mu, r = r, sigma = sigma,
  nsteps = nsteps, sizesteps = dt, npaths = npaths,
  theseed = theseed, isPmeasure = FALSE
)

StockPathsQantithetic <- SimulateStockPathsAntithetic(
  S = S0, mu = mu, r = r, sigma = sigma,
  nsteps = nsteps, sizesteps = dt, npaths = npaths,
  theseed = theseed, isPmeasure = FALSE
)

# d)

# Extraire les dernières valeurs ST
ST_P <- StockPathsP[, ncol(StockPathsP)]
ST_P_anti <- StockPathsPantithetic[, ncol(StockPathsPantithetic)]
ST_Q <- StockPathsQ[, ncol(StockPathsQ)]
ST_Q_anti <- StockPathsQantithetic[, ncol(StockPathsQantithetic)]

# Fonction pour calculer IC
confidence_interval_95 <- function(x) {
  m <- mean(x)
  s <- sd(x)
  n <- length(x)
  error_margin <- 1.96 * s / sqrt(n)
  c(m - error_margin, m + error_margin)
}

# IC pour E_P[ST] avec trajectoires normales
IC_P <- confidence_interval_95(ST_P)
# IC pour E_P[ST] avec trajectoires antithétiques
IC_P_anti <- confidence_interval_95(ST_P_anti)

# IC pour E_Q[ST] avec trajectoires normales
IC_Q <- confidence_interval_95(ST_Q)
# IC pour E_Q[ST] avec trajectoires antithétiques
IC_Q_anti <- confidence_interval_95(ST_Q_anti)

# Affichage des résultats
cat("Estimation E_P[ST] (classique) :", mean(ST_P), "\nIC 95% :", IC_P, "\n")
cat("Estimation E_P[ST] (antithétique) :", mean(ST_P_anti), "\nIC 95% :", IC_P_anti, "\n")

cat("Estimation E_Q[ST] (classique) :", mean(ST_Q), "\nIC 95% :", IC_Q, "\n")
cat("Estimation E_Q[ST] (antithétique) :", mean(ST_Q_anti), "\nIC 95% :", IC_Q_anti, "\n")

# Calcul des valeurs exactes
exact_Ep <- S0 * exp(mu * T)
exact_Eq <- S0 * exp(r * T)
cat("Valeur exacte E_P[ST] :", exact_Ep, "\n")
cat("Valeur exacte E_Q[ST] :", exact_Eq, "\n")

#Exercice 5

# a)

# Supposez que les matrices StockPathsQ et StockPathsQantithetic sont déjà en mémoire.
# Dimensions : M x 27, M = 10 000 par exemple.
# StockPathsQ[i, ] contient la trajectoire i, des temps t0 à t26.

# Paramètres
r <- 0.02
T <- 0.5
K <- 100

# Nombre de trajectoires
M <- nrow(StockPathsQ)

# Calcul de la moyenne arithmétique pour chaque trajectoire
# rowMeans calcule la moyenne de chaque ligne
A_T_Q <- rowMeans(StockPathsQ)
A_T_Q_anti <- rowMeans(StockPathsQantithetic)

# Calcul du payoff
payoff_Q <- pmax(A_T_Q - K, 0)
payoff_Q_anti <- pmax(A_T_Q_anti - K, 0)

# Actualisation des payoffs
discount_factor <- exp(-r * T)
discounted_payoff_Q <- discount_factor * payoff_Q
discounted_payoff_Q_anti <- discount_factor * payoff_Q_anti

# Estimation du prix par Monte Carlo
price_estimate_Q <- mean(discounted_payoff_Q)
price_estimate_Q_anti <- mean(discounted_payoff_Q_anti)

# Écart-type empirique
std_Q <- sd(discounted_payoff_Q)
std_Q_anti <- sd(discounted_payoff_Q_anti)

# Intervalle de confiance 95%
# IC : estimate ± 1.96 * (std / sqrt(M))
error_margin_Q <- 1.96 * std_Q / sqrt(M)
error_margin_Q_anti <- 1.96 * std_Q_anti / sqrt(M)

IC_Q <- c(price_estimate_Q - error_margin_Q, price_estimate_Q + error_margin_Q)
IC_Q_anti <- c(price_estimate_Q_anti - error_margin_Q_anti, price_estimate_Q_anti + error_margin_Q_anti)

# Affichage des résultats
cat("Estimation du prix de l'option asiatique sans antithétique :\n")
cat("Prix estimé :", price_estimate_Q, "\n")
cat("IC 95% :", IC_Q, "\n\n")

cat("Estimation du prix de l'option asiatique avec antithétique :\n")
cat("Prix estimé :", price_estimate_Q_anti, "\n")
cat("IC 95% :", IC_Q_anti, "\n")

# c)

# Suppose que StockPathsQ et StockPathsQantithetic sont déjà disponibles.
# Paramètres
r <- 0.02
T <- 0.5
K <- 100
dt <- 1/52
n <- 26 # Nombre de pas (hors temps initial) => total de 27 points: t0 à t26
M <- nrow(StockPathsQ) # nombre de trajectoires
S0 <- 100 # Comme dans l'exercice 4
mu <- 0.07  # juste pour rappel, pas utilisé dans la mesure Q
sigma <- 0.20

# 1. Calcul de A_T^{(ari)} pour chaque trajectoire sous Q
A_T_Q <- rowMeans(StockPathsQ)
A_T_Q_anti <- rowMeans(StockPathsQantithetic)

# 2. Payoff H^{(ari)} = max(A_T - K, 0)
payoff_Q <- pmax(A_T_Q - K, 0)
payoff_Q_anti <- pmax(A_T_Q_anti - K, 0)

# 3. Actualisation
discount_factor <- exp(-r * T)
X_Q <- discount_factor * payoff_Q
X_Q_anti <- discount_factor * payoff_Q_anti

# 4. Calcul de E_Q[A_T^{(ari)}]:
# Formellement : E_Q[A_T^{(ari)}] = (S0/(n)) * e^{r dt}*(1 - e^{r n dt})/(1 - e^{r dt})
# On sait qu'il y a n=26 pas, dt=1/52, donc n dt = T = 0.5.
EQ_A <- (S0 / n) * (exp(r*dt) * (1 - exp(r*n*dt)) / (1 - exp(r*dt)))

# 5. Calcul du c^*
# c^* = - Cov(X, Y) / Var(Y), où X = e^{-rT}H^{(ari)} et Y = A_T^{(ari)}
Y <- A_T_Q
X <- X_Q

cov_XY <- cov(X, Y)
var_Y <- var(Y)
c_star <- - cov_XY / var_Y

# 6. Estimateur par variable de contrôle
# Estimateur = mean( X + c_star * (Y - E_Q[A_T^{(ari)}]) )
X_cv <- X + c_star * (Y - EQ_A)
price_cv <- mean(X_cv)

# IC 95%
std_cv <- sd(X_cv)
error_margin_cv <- 1.96 * std_cv / sqrt(M)
IC_cv <- c(price_cv - error_margin_cv, price_cv + error_margin_cv)

# Même chose pour les trajectoires antithétiques :
Y_anti <- A_T_Q_anti
X_anti <- X_Q_anti

cov_XY_anti <- cov(X_anti, Y_anti)
var_Y_anti <- var(Y_anti)
c_star_anti <- - cov_XY_anti / var_Y_anti

X_cv_anti <- X_anti + c_star_anti * (Y_anti - EQ_A)
price_cv_anti <- mean(X_cv_anti)

std_cv_anti <- sd(X_cv_anti)
error_margin_cv_anti <- 1.96 * std_cv_anti / sqrt(M)
IC_cv_anti <- c(price_cv_anti - error_margin_cv_anti, price_cv_anti + error_margin_cv_anti)

# Affichage des résultats
cat("Résultats par variable de contrôle (sans antithétique):\n")
cat("c* =", c_star, "\n")
cat("Prix estimé:", price_cv, "\n")
cat("IC 95%:", IC_cv, "\n\n")

cat("Résultats par variable de contrôle (avec antithétique):\n")
cat("c* =", c_star_anti, "\n")
cat("Prix estimé:", price_cv_anti, "\n")
cat("IC 95%:", IC_cv_anti, "\n")