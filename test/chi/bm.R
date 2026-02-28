library(ape)
library(PCMBase)

# ultrametric tree of 3 tips:
newick <- "((K:0.08,C:0.08):0.05,H:0.13);"
# newick <- "(K:0.08,C:0.08,H:0.08):0;"
tree <- PCMTree(read.tree(text = newick))

# trait data:
X <- matrix(c(0.1,  NaN,   NA,
             -0.2,  NaN,  1.2,
               NA,  0.4,  0.7), nrow = 3)
colnames(X) <- c("K", "C", "H")
# correlation matrix:
R <- matrix(c(1.0,  0.3,  0.0,
              0.3,  1.0, -0.5,
              0.0, -0.5,  1.0), nrow = 3)

model.BM <- PCM(model = "BM", k = nrow(X))
model.BM$X0[] <- c(NA, NA, NA)
model.BM$Sigma_x[,,1] <- t(chol(R))
# note: the chol function gives an upper triangular matrix,
# i.e., L', and R = L'L

# Likelihood:
PCMLik(X, tree, model.BM)

options(digits = 3)
print(PCMTable(model.BM, removeUntransformed = FALSE),
      xtable = TRUE, type = "latex")

# tracing the likelihood calculation
porder <- c("K", "C", "H", "5", "4")
traceTable <- PCMLikTrace(X, tree, model.BM)
traceTable <- traceTable[match(porder, traceTable$i)]

cat(FormatTableAsLatex(traceTable[, list(i, t_i, k_i, X_i, V_i, Phi_i, A_i)]))
cat(FormatTableAsLatex(traceTable[, list(C_i, E_i, L_i, m_i, r_i)]))
cat(FormatTableAsLatex(traceTable[, list(j, i, `L_{ji}`, `m_{ji}`, `r_{ji}`,
                                         `\\hat{X}_i`, `\\ell\\ell_i`)]))
L0 <- traceTable[5,]$L_i[[1]]
m0 <- traceTable[5,]$m_i[[1]]
r0 <- traceTable[5,]$r_i[[1]]
x0 <- -0.5 * solve(L0) %*% m0
ll <- t(x0) %*% L0 %*% x0 + t(x0) %*% m0 + r0


# original data:
M <- matrix(c(-1.3595, -0.8086,
              -0.9219, -0.8869,
              -0.1170, -1.0917), nrow = 3, byrow = TRUE)
R <- matrix(c(1, 0.25, 0.25, 1), nrow = 2)
logdetR <- log(det(R))
# transformed data:
Z <-  M %*% solve(chol(R))
