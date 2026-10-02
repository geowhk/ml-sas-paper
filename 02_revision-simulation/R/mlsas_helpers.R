as_dgRMatrix_listw <- function(listw) {
  neighbours <- listw$neighbours
  weights <- listw$weights
  n <- length(neighbours)

  rows <- rep(seq_len(n), lengths(neighbours))
  cols <- unlist(neighbours, use.names = FALSE)
  vals <- unlist(weights, use.names = FALSE)
  keep <- cols > 0

  Matrix::sparseMatrix(
    i = rows[keep],
    j = cols[keep],
    x = vals[keep],
    dims = c(n, n)
  )
}

compute_path_proximity <- function(
    data,
    n_tree = 500,
    workers = 3,
    validate = FALSE,
    mtry = NULL,
    nodesize = NULL
) {
  
  # 입력 자료는 data.frame 형태여야 한다.
  n <- nrow(data)
  
  # ---- 1. 랜덤 포레스트 모형 적합(Addcl1 방식) ----
  synthetic <- as.data.frame(lapply(data, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", nrow(data)), rep("synthetic", nrow(synthetic))))
  
  rf_args <- list(
    x = data_all,
    y = y_all,
    ntree = n_tree
  )
  if (!is.null(mtry)) {
    rf_args$mtry <- as.integer(mtry)
  }
  if (!is.null(nodesize)) {
    rf_args$nodesize <- as.integer(nodesize)
  }
  rf <- do.call(randomForest, rf_args)
  
  # ---- 2. 각 나무에서 관측치별 분기 경로를 계산하는 함수 ----
  get_paths <- function(tree_num) {
    tree <- getTree(rf, k = tree_num, labelVar = FALSE)
    paths <- vector("list", 2 * n)
    
    for (i in 1:(2 * n)) {
      path <- c()
      node <- 1
      repeat {
        path <- c(path, node)
        split_var <- tree[node, "split var"]
        if (split_var == 0) break
        
        split_val <- tree[node, "split point"]
        split_var_name <- colnames(data_all)[split_var]
        value <- data_all[i, split_var_name]
        
        split_column <- data_all[[split_var_name]]

        if (is.factor(split_column) && !is.ordered(split_column)) {
          split_bits <- intToBits(split_val)
          value_level_index <- as.integer(value)
          
          if (split_bits[value_level_index] == as.raw(1)) {
            node <- tree[node, "left daughter"]
          } else {
            node <- tree[node, "right daughter"]
          }
          
        } else {
          split_value <- if (is.ordered(split_column)) as.integer(value) else value
          if (split_value <= split_val) {
            node <- tree[node, "left daughter"]
          } else {
            node <- tree[node, "right daughter"]
          }
        }
      }
      paths[[i]] <- path
    }
    return(paths)
  }
  
  # ---- [검증 코드] ----
  if (validate == TRUE) {
    message("Validation step is running...")
    predicted_nodes <- attr(predict(rf, data_all, nodes = TRUE), "nodes")
    all_paths <- future_lapply(1:n_tree, get_paths, future.packages = "randomForest")
    
    mismatch_found <- FALSE
    for (t in 1:n_tree) {
      truth_nodes_for_tree_t <- predicted_nodes[, t]
      user_paths_for_tree_t <- all_paths[[t]]
      user_final_nodes_for_tree_t <- sapply(user_paths_for_tree_t, function(p) p[length(p)])
      if (!all(truth_nodes_for_tree_t == user_final_nodes_for_tree_t)) {
        mismatch_found <- TRUE
        break
      }
    }
    if (mismatch_found) {
      warning("Validation Failed!")
    } else {
      message("Validation Success: 'get_paths' function works.")
    }
  } else {
    all_paths <- future_lapply(1:n_tree, get_paths, future.packages = "randomForest")
  }
  
  # ---- [검증 종료] ----
  
  
  # ---- 4. 근접도 행렬 병렬 계산 ----
  message("Calculating path proximity matrix in parallel...")
  proximity <- matrix(0, nrow = n, ncol = n)
  
  tree_jaccard <- function(b) {
    paths <- all_paths[[b]]
    mtx <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      Ti <- paths[[i]]
      for (j in i:n) {
        Tj <- paths[[j]]
        intersect_len <- length(intersect(Ti, Tj))
        union_len <- length(union(Ti, Tj))
        jaccard <- intersect_len / union_len
        mtx[i, j] <- jaccard
        if (i != j) mtx[j, i] <- jaccard
      }
    }
    return(mtx)
  }
  
  
  mtx_list <- future_lapply(1:n_tree, tree_jaccard, future.globals = c("all_paths", "n"))
  
  
  # ---- 5. 나무별 결과를 합산하고 정규화 ----
  for (m in mtx_list) {
    proximity <- proximity + m
  }
  proximity <- proximity / n_tree
  
  rownames(proximity) <- rownames(data)
  colnames(proximity) <- rownames(data)
  
  return(proximity)
}

permute_cols_rowwise <- function(Q) {
  n <- nrow(Q)
  Qp <- Q
  for (i in seq_len(n)) {
    idx <- setdiff(seq_len(n), i)               # 대각 제외
    Qp[i, idx] <- Q[i, sample(idx)]             # 행 i에서 열만 섞기
  }
  diag(Qp) <- 0
  Qp
}

permute_mantel <- function(Q) {
  n <- nrow(Q)
  pi <- sample.int(n)       # 1:n의 무작위 순열
  Qp <- Q[pi, pi]           # 행과 열 모두 동일한 순서로 섞기
  diag(Qp) <- 0             # 자기 자신과의 관계는 0으로 (선택)
  Qp
}

global.mlsas <- function(nb, S, style="W", diag.zero=TRUE) {
  style <- match.arg(style, c("B","C","S","W","U","V"))
  mf <- match.call(expand.dots = FALSE)
  
  # 공간 가중치 행렬의 대각 원소를 0이 아닌 값으로 둘지 처리한다.
  if (style=="V") diag.zero <- FALSE
  if (diag.zero == FALSE) nb <- nb.self(nb)
  if (style=="V") {listw <- nb2listv(nb,beta=0.5)}
  else {listw <- nb2listw(nb,style=style)}
  
  V <- as.matrix(as_dgRMatrix_listw(listw))

  n <- nrow(S)
  
  diag(V) <- 0
  diag(S) <- 0
  P <- V/sum(V)
  Q <- S/sum(S)

  mlsas <- n*(n-1)*sum(P*Q)
  return(mlsas)
}

test.global.mlsas <- function(nb, S, style="W", method="randomization", alternative="two.sided", diag.zero=TRUE, nperm = 999){
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  style <- match.arg(style, c("B","C","S","W","U","V"))
  method <- match.arg(method, c("randomization", "permutation"))
  mf <- match.call(expand.dots = FALSE)
  
  # 공간 가중치 행렬의 대각 원소를 0이 아닌 값으로 둘지 처리한다.
  if (style=="V") diag.zero <- FALSE
  if (diag.zero == FALSE) nb <- nb.self(nb)
  if (style=="V") {listw <- nb2listv(nb,beta=0.5)}
  else {listw <- nb2listw(nb,style=style)}
  # 필요하면 listw2sn(listw)로 변환해 확인할 수 있다.
  
  V <- as.matrix(as_dgRMatrix_listw(listw))
  n <- nrow(S)
  diag(V) <- 0
  diag(S) <- 0
  P <- n*(n-1)*V/sum(V)
  Q <- S/sum(S)

  mlsas <- sum(P*Q)
  
  if(method == "randomization"){
  
    P <- 0.5*(P+t(P))
    Q <- 0.5*(Q+t(Q))
    P.off <- P
    diag(P.off) <- 0
    Q.off <- Q
    diag(Q.off) <- 0
    P.on <- diag(P)
    Q.on <- diag(Q)
    
    F0.off <- sum(P.off)
    F1.off <- sum(P.off^2)
    F2.off <- sum(rowSums(P.off)^2)
    F0.on <- sum(P.on)
    F1.on <- sum(P.on^2)
    F0.all <- sum(P)
    F1.all <- sum(P^2)
    F2.all <- sum(rowSums(P)^2)
    F2.on <- sum(P.on*rowSums(P))
    G0.off <- sum(Q.off)
    G1.off <- sum(Q.off^2)
    G2.off <- sum(rowSums(Q.off)^2)
    G0.on <- sum(Q.on)
    G1.on <- sum(Q.on^2)
    G2.on <- sum(Q.on*rowSums(Q))
    G0.all <- sum(Q)
    G1.all <- sum(Q^2)
    G2.all <- sum(rowSums(Q)^2)
    
    E.I.off <- (F0.off*G0.off)/(n*(n-1))
    E.I.on <- (F0.on*G0.on)/n
    E.I <- E.I.off + E.I.on
    
    VAR.I.off <- ((2*F1.off*G1.off)/(n*(n-1)))+((4*(F2.off-F1.off)*(G2.off-G1.off))/(n*(n-1)*(n-2)))+ ((((F0.off^2)+(2*F1.off)-(4*F2.off))*((G0.off^2)+(2*G1.off)-(4*G2.off)))/(n*(n-1)*(n-2)*(n-3)))-(E.I.off^2)
    VAR.I.on <- (F1.on*G1.on/n)+(((F0.on^2-F1.on)*(G0.on^2-G1.on))/(n*(n-1)))-(E.I.on^2)
    COV.1 <- ((F2.all-F1.on-F2.off)*(G2.all-G1.on-G2.off))/(2*n*(n-1))
    COV.2 <- (((F0.on*F0.off)-(F2.all-F1.on-F2.off))*((G0.on*G0.off)-(G2.all-G1.on-G2.off)))/(n*(n-1)*(n-2))
    COV.I <- (COV.1+COV.2)-(E.I.off*E.I.on)
    COV.I.2 <- 2*COV.I
    VAR.I <- VAR.I.off + VAR.I.on + COV.I.2
    SD.I <- sqrt(VAR.I)
    Z.I <- (mlsas - E.I)/SD.I

    if (alternative == "two.sided"){
      P.I <- 2 * pnorm(abs(Z.I), lower.tail = FALSE)}
    else if (alternative == "greater"){
      P.I <- pnorm(Z.I, lower.tail = FALSE)}
    else {P.I <- pnorm(Z.I)}
    
    statistic <- Z.I
    attr(statistic, "names") <- "ML-SAS statistic standard deviate"
    p.value <- P.I
    estimate <- c(mlsas, E.I, VAR.I)
    attr(estimate, "names") <- c("Observed ML-SAS", "Expectation","Variance")
    method <- "Global ML-SAS test under randomization assumption"
    
    if (diag.zero==FALSE) style <- paste(style, "*",sep="")
    data.name <- paste(deparse(mf[[2]]), "\nweights: ",  deparse(mf[[3]]), ", style: ", style, "\n", sep = "")
    res <- list(statistic = statistic, p.value = p.value, estimate = estimate,
                method = method, alternative = alternative, data.name = data.name)
    class(res) <- "htest"
    res
  }
  else if(method == "permutation"){
    perm_stats <- numeric(nperm)
    for (b in seq_len(nperm)) {
      Sp <- permute_mantel(S)
      perm_stats[b] <- global.mlsas(nb=nb, Sp, style=style, diag.zero=diag.zero)
    }
    E.I <- mean(perm_stats)
    VAR.I <- var(perm_stats)
    SKEW.I <- (sum((perm_stats - E.I)^3)/nperm)/((sum((perm_stats - E.I)^2)/nperm)^(3/2))
    KURT.I <- (sum((perm_stats - E.I)^4)/nperm)/((sum((perm_stats - E.I)^2)/nperm)^2) - 3
    SD.I <- sqrt(VAR.I)
    Z.I <- (mlsas - E.I)/SD.I
    
    if (alternative == "two.sided"){
      P.I <- 2 * pnorm(abs(Z.I), lower.tail = FALSE)}
    else if (alternative == "greater"){
      P.I <- pnorm(Z.I, lower.tail = FALSE)}
    else {P.I <- pnorm(Z.I)}
    
    p_perm_high <- (sum(perm_stats >= mlsas) + 1) / (nperm + 1)
    p_perm_low <- (sum(perm_stats < mlsas) + 1) / (nperm + 1)
    p_perm_two <- pmin(p_perm_high, p_perm_low)
    
    if (alternative == "two.sided"){
      p_perm <- p_perm_two * 2}
    else if (alternative == "greater"){
      p_perm <- p_perm_high}
    else {p_perm <- p_perm_low}
    
    statistic <- Z.I
    attr(statistic, "names") <- "ML-SAS statistic standard deviate"
    p.value <- P.I
    estimate <- c(mlsas, E.I, VAR.I, SKEW.I, KURT.I)
    attr(estimate, "names") <- c("Observed ML-SAS", "Expectation","Variance", "Skewness", "Kurtosis")
    method <- "Global ML-SAS Permutation test"
    
    if (diag.zero==FALSE) style <- paste(style, "*",sep="")
    data.name <- paste(deparse(mf[[2]]), "\nweights: ",  deparse(mf[[3]]), ", style: ", style, "\n", sep = "")
    res <- list(statistic = statistic, p.value = p.value, estimate = estimate,
                method = method, alternative = alternative, data.name = data.name)
    class(res) <- "htest"
    res     
  }
}

local.mlsas <- function(id, nb, S, style="W", diag.zero=TRUE){

  style <- match.arg(style, c("B","C","S","W","U","V"))
  mf <- match.call(expand.dots = FALSE)
  
  n <- nrow(S)
  
  if (style=="V") diag.zero <- FALSE
  if (diag.zero == FALSE) nb <- nb.self(nb)
  if (style=="V") {listw <- nb2listv(nb,beta=0.5)}
  else {listw <- nb2listw(nb,style=style)}
  
  V <- as.matrix(as_dgRMatrix_listw(listw))
  diag(V) <- 0
  diag(S) <- 0
  P <- n*(n-1)*V/sum(V)
  Q <- S/sum(S)

  mlsas.vec <- n*rowSums(P*Q)
  
  res.local.mlsas <- data.frame(id=id, local.mlsas=mlsas.vec)
  return(res.local.mlsas)
}

test.local.mlsas <- function(id, nb, S, style="W", sig.levels=c(0.05), method="total", alternative="two.sided", diag.zero=TRUE, p.adj = "fdr", nperm = 999){

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  style <- match.arg(style, c("B","C","S","W","U","V"))
  method <- match.arg(method, c("conditional", "total", "permute_cond", "permute_total"))
  p.adj <- match.arg(p.adj, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  mf <- match.call(expand.dots = FALSE)
  
  n <- nrow(S)

  if (style=="V") diag.zero <- FALSE
  if (diag.zero == FALSE) nb <- nb.self(nb)
  if (style=="V") {listw <- nb2listv(nb,beta=0.5)}
  else {listw <- nb2listw(nb,style=style)}
  
  V <- as.matrix(as_dgRMatrix_listw(listw))
  diag(V) <- 0
  diag(S) <- 0
  P <- n*(n-1)*V/sum(V)
  Q <- S/sum(S)
  
  local.moran.vec <- n*rowSums(P*Q)

  # 결과 저장용 객체

  E.I.vec <- vector()
  VAR.I.vec <- vector()

  # 각 관측치별 계산
  
  Q <- 0.5*(Q+t(Q))

  if (method=="conditional" | method=="total"){
    
  for (i in 1:n){
    vi.vec <- listw$weights[[i]]
    vi.s <- sum(vi.vec)
    # 기대값과 분산 계산
    if (method=="conditional") {
      p_i <- P[i, ]
      p_i <- p_i[-i]*n
      q_i <- Q[i, ]
      q_i <- q_i[-i]
      E.I <- (1/(n-1))*sum(p_i)*sum(q_i)
      VAR.I <- (1/(n-2))*sum((p_i - mean(p_i))^2)*sum((q_i - mean(q_i))^2)
    }
    else if (method=="total") {
      P_i <- P
      P_i[-i, ] <- 0
      P_i <- P_i * n
      P_i <- 0.5*(P_i+t(P_i))
      P.off <- P_i
      diag(P.off) <- 0
      Q.off <- Q
      diag(Q.off) <- 0
      P.on <- diag(P_i)
      Q.on <- diag(Q)
      
      F0.off <- sum(P.off)
      F1.off <- sum(P.off^2)
      F2.off <- sum(rowSums(P.off)^2)
      F0.on <- sum(P.on)
      F1.on <- sum(P.on^2)
      F0.all <- sum(P)
      F1.all <- sum(P^2)
      F2.all <- sum(rowSums(P)^2)
      F2.on <- sum(P.on*rowSums(P))
      G0.off <- sum(Q.off)
      G1.off <- sum(Q.off^2)
      G2.off <- sum(rowSums(Q.off)^2)
      G0.on <- sum(Q.on)
      G1.on <- sum(Q.on^2)
      G2.on <- sum(Q.on*rowSums(Q))
      G0.all <- sum(Q)
      G1.all <- sum(Q^2)
      G2.all <- sum(rowSums(Q)^2)
      
      E.I.off <- (F0.off*G0.off)/(n*(n-1))
      E.I.on <- (F0.on*G0.on)/n
      E.I <- E.I.off + E.I.on
      
      VAR.I.off <- ((2*F1.off*G1.off)/(n*(n-1)))+((4*(F2.off-F1.off)*(G2.off-G1.off))/(n*(n-1)*(n-2)))+ ((((F0.off^2)+(2*F1.off)-(4*F2.off))*((G0.off^2)+(2*G1.off)-(4*G2.off)))/(n*(n-1)*(n-2)*(n-3)))-(E.I.off^2)
      VAR.I.on <- (F1.on*G1.on/n)+(((F0.on^2-F1.on)*(G0.on^2-G1.on))/(n*(n-1)))-(E.I.on^2)
      COV.1 <- ((F2.all-F1.on-F2.off)*(G2.all-G1.on-G2.off))/(2*n*(n-1))
      COV.2 <- (((F0.on*F0.off)-(F2.all-F1.on-F2.off))*((G0.on*G0.off)-(G2.all-G1.on-G2.off)))/(n*(n-1)*(n-2))
      COV.I <- (COV.1+COV.2)-(E.I.off*E.I.on)
      COV.I.2 <- 2*COV.I
      VAR.I <- VAR.I.off + VAR.I.on + COV.I.2
    }
    E.I.vec[i] <- E.I
    VAR.I.vec[i] <- VAR.I

  }
  SD.I.vec <- sqrt(VAR.I.vec)
  Z.I.vec <- (local.moran.vec-E.I.vec)/SD.I.vec
  
  # p-value 계산
  # P.L.vec <- numeric(n)
  if (alternative=="two.sided") {
    P.I.vec <- sapply(Z.I.vec, function(x) {2 * pnorm(abs(x), lower.tail = FALSE)})}
  else if (alternative=="greater") {
    P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = FALSE)})}
  else {
    P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = TRUE)})}
  
  
    local.Moran <- data.frame(id=id, local.Moran=local.moran.vec, E=E.I.vec, VAR=VAR.I.vec,
                              SD=SD.I.vec, Z=Z.I.vec, P=P.I.vec, 
                              P.adj = p.adjust(P.I.vec, method = p.adj))
    return(local.Moran)
  }
  if(method=="permute_total"){
    perm_mat <- matrix(NA_real_, nrow = n, ncol = nperm)
    for (b in seq_len(nperm)) {
      Sp <- permute_mantel(S)
      perm_mat[, b] <- local.mlsas(1:n, nb=nb, Sp, style=style, diag.zero=diag.zero)[, 2]
    }
    E.I.vec <- rowMeans(perm_mat)
    VAR.I.vec <- apply(perm_mat, 1, var)
    SD.I.vec <- sqrt(VAR.I.vec)
    Z.I.vec <- (local.moran.vec-E.I.vec)/SD.I.vec
    if (alternative=="two.sided") {
      P.I.vec <- sapply(Z.I.vec, function(x) {2 * pnorm(abs(x), lower.tail = FALSE)})}
    else if (alternative=="greater") {
      P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = FALSE)})}
    else {
      P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = TRUE)})}

    p_perm_high <- (rowSums(perm_mat >= local.moran.vec) + 1) / (nperm + 1)
    p_perm_low <- (rowSums(perm_mat < local.moran.vec) + 1) / (nperm + 1)
    p_perm_two <- pmin(p_perm_high, p_perm_low)
    
    if (alternative == "two.sided"){
      p_perm <- p_perm_two * 2}
    else if (alternative == "greater"){
      p_perm <- p_perm_high}
    else {p_perm <- p_perm_low}
            
    # greater 대립가설 방향의 permutation p-value 계산식

    local.Moran <- data.frame(id=id, local.Moran=local.moran.vec, E=E.I.vec, VAR=VAR.I.vec,
                              SD=SD.I.vec, Z=Z.I.vec, P=P.I.vec, p_perm = p_perm, 
                              P.adj = p.adjust(P.I.vec, method = p.adj), 
                              p_perm_adj = p.adjust(p_perm, method = p.adj))
    return(local.Moran)
  }
  if(method=="permute_cond"){
    perm_mat <- matrix(NA_real_, nrow = n, ncol = nperm)
    for (b in seq_len(nperm)) {
      Sp <- permute_cols_rowwise(S)
      perm_mat[, b] <- local.mlsas(1:n, nb=nb, Sp, style=style, diag.zero=diag.zero)[, 2]
    }
    E.I.vec <- rowMeans(perm_mat)
    VAR.I.vec <- apply(perm_mat, 1, var)
    SD.I.vec <- sqrt(VAR.I.vec)
    Z.I.vec <- (local.moran.vec-E.I.vec)/SD.I.vec
    if (alternative=="two.sided") {
      P.I.vec <- sapply(Z.I.vec, function(x) {2 * pnorm(abs(x), lower.tail = FALSE)})}
    else if (alternative=="greater") {
      P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = FALSE)})}
    else {
      P.I.vec <- sapply(Z.I.vec, function(x) {pnorm(x, lower.tail = TRUE)})}
    
    p_perm_high <- (rowSums(perm_mat >= local.moran.vec) + 1) / (nperm + 1)
    p_perm_low <- (rowSums(perm_mat < local.moran.vec) + 1) / (nperm + 1)
    p_perm_two <- pmin(p_perm_high, p_perm_low)
    
    if (alternative == "two.sided"){
      p_perm <- p_perm_two * 2}
    else if (alternative == "greater"){
      p_perm <- p_perm_high}
    else {p_perm <- p_perm_low}
    
    # greater 대립가설 방향의 permutation p-value 계산식
    
    local.Moran <- data.frame(id=id, local.Moran=local.moran.vec, E=E.I.vec, VAR=VAR.I.vec,
                              SD=SD.I.vec, Z=Z.I.vec, 
                              P=P.I.vec, P.adj = p.adjust(P.I.vec, method = p.adj), 
                              p_perm = p_perm, p_perm_adj = p.adjust(p_perm, method = p.adj))
    return(local.Moran)
  }
}

snu_local_geary <- function(localgeary_res, p.adj.type = "fdr", q.value = 0.05){
  p.adj.type <- match.arg(p.adj.type, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  res_geary <- tibble(
    ci = data.frame(local_geary)[, 1],
    E.ci = data.frame(attr(local_geary, "pseudo-p"))[, 1],
    Var.ci = data.frame(attr(local_geary, "pseudo-p"))[, 2],
    Z.ci = data.frame(attr(local_geary, "pseudo-p"))[, 3],
    P.ci =  data.frame(attr(local_geary, "pseudo-p"))[, 4],
    P.ci.sim =  data.frame(attr(local_geary, "pseudo-p"))[, 5],
    P.ci.sim.fold =  data.frame(attr(local_geary, "pseudo-p"))[, 6],
    class = data.frame(attr(local_geary, "cluster"))[, 1]
  ) |> 
    mutate(
      p.adj = p.adjust(P.ci, method = p.adj.type), # p.ci는 cluster detection을 못함.
      p.adj.sim = p.adjust(P.ci.sim, method = p.adj.type), 
      class.sig = if_else(
        p.adj > q.value, "Not Significant", class
      ),
      class.sig.sim = if_else(
        p.adj.sim > q.value, "Not Significant", class
      ),
      class.sig = fct(class.sig, levels = c("High-High", "Low-Low", "Other Positive", "Negative", "Not Significant")),
      class.sig.sim = fct(
        class.sig.sim, levels = c("High-High", "Low-Low", "Other Positive", "Negative", "Not Significant")
      )
    )
  return(res_geary)  
}
