train.corels <- function(tdata,pos_sign="1", neg_sign="0",rule_minlen=1,rule_maxlen=1,minsupport_pos=0.10,minsupport_neg=0.10,curiosity_policy=2,max_nodes=10000, regularization=0.01, verbosity="progress", map_type=1, ablation=0, calculate_size=0, latex_out=1) {
    train_from_file <- function(data_fname,label_fname,meta_fname="",curiosity_policy=2,max_nodes=10000, regularization=0.01, verbosity="progress", map_type=1, ablation=0, calculate_size=0, latex_out=0) {
        # ugly hack to get Rcpp to be able to recognize the parameters
        params_list <- c(toString(curiosity_policy), format(max_nodes, scientific=FALSE), toString(regularization), toString(verbosity), toString(map_type), toString(ablation), toString(calculate_size), toString(latex_out))
        rs <- .Call('corels_train', PACKAGE='corels', params_list, data_fname, label_fname, meta_fname)
        rs
    }

    pos_data <- tdata[tdata$label==pos_sign,]
    neg_data <- tdata[tdata$label==neg_sign,]

    pos_rules <- arules::eclat(subset(pos_data, select =- label), parameter = list(minlen=rule_minlen, maxlen=rule_maxlen, supp=minsupport_pos))
    neg_rules <- arules::eclat(subset(neg_data, select =- label), parameter = list(minlen=rule_minlen, maxlen=rule_maxlen, supp=minsupport_neg))

    pos_featurenames <- attributes(attributes(pos_rules)$items)$itemInfo$labels
    pos_rulenames <- as(pos_rules, "data.frame")$items #inspect(pos_rules)$items
    pos_mat <- attributes(attributes(pos_rules)$items)$data
    neg_featurenames <- attributes(attributes(neg_rules)$items)$itemInfo$labels
    neg_rulenames <- as(neg_rules, "data.frame")$items #inspect(neg_rules)$items
    neg_mat <- attributes(attributes(neg_rules)$items)$data

    pos_mat2 <- rbind(pos_mat, matrix(0, length(setdiff(neg_featurenames, pos_featurenames)), ncol(pos_mat)))
    neg_mat2 <- rbind(neg_mat, matrix(0, length(setdiff(pos_featurenames, neg_featurenames)), ncol(neg_mat)))

    pos_combined_featurenames <- c(as.character(pos_featurenames),as.character(setdiff(neg_featurenames, pos_featurenames)))
    neg_combined_featurenames <- c(as.character(neg_featurenames),as.character(setdiff(pos_featurenames, neg_featurenames)))
    featurenames <- sort(pos_combined_featurenames, index.return=TRUE) # all the features
    # indices for the positive features
    pos_idx <- featurenames$ix
    featurenames <- featurenames$x
    # indices for the negative features
    neg_idx <- sort(neg_combined_featurenames, index.return=TRUE)$ix

    pos_mat3 <- pos_mat2[pos_idx, ]
    neg_mat3 <- neg_mat2[neg_idx, ]
    # now we had the rows correct. let's fix the columns (rules)
    rulenames <- c(as.character(pos_rulenames), as.character(neg_rulenames))
    idx <- order(rulenames)[!duplicated(sort(rulenames))]
    rulenames <- rulenames[idx]
    # get the columns correct for feature_rule matrix
    mat <- as.matrix(cbind(pos_mat3, neg_mat3)[, idx])

    mat_data_feature <- get_data_feature_mat(tdata, featurenames)
    # get the data_rule matrix by multiplying data_feature and feature_rule matrices
    mat_data_rules <- mat_data_feature %*% mat
    mat_data_rules <- t(t(mat_data_rules)>=c(colSums(mat)))+0

    write.table(as.matrix(t(mat_data_rules)), file='tdata_R.out', sep=' ', row.names=rulenames, col.names=FALSE, quote=FALSE)
    label <- t(cbind((tdata$label==neg_sign) +0, (tdata$label==pos_sign) +0))
    write.table(as.matrix(label), file='tdata_R.label', sep=' ', row.names=c("{label=0}", "{label=1}"), col.names=FALSE, quote=FALSE)

    label <- t(label)
    r_samples <- character(dim(mat_data_rules)[1])
    for (i in 1:dim(mat_data_rules)[1]) {
        r_samples[i] <- paste(mat_data_rules[i,], collapse='')
    }

    df = data.frame(r_samples, label[,1], list(rep(1,dim(mat_data_rules)[1])))
    colnames(df) <- c("obs", "label", "count")
    agg <- aggregate(cbind(label=df$label, count=df$count), by=list(obs=df$obs), FUN=sum)
    ord <- agg[order(-agg$count),]
    ml <- as.integer(ord$label < ord$count-ord$label)
    d <- data.frame(obs=agg$obs, ml=ml)

    ind <- integer(dim(df)[1])
    for (i in 1:dim(df)[1]) {
        ind[i] <- as.integer(df$label[i] == d[which(d$obs==df$obs[i]),]$ml)
    }
    write.table(t(as.matrix(ind)), file='tdata_R.minor', sep=' ', row.names=c("{group_minority}"), col.names=FALSE, quote=FALSE)

    rs <- train_from_file('tdata_R.out', 'tdata_R.label', meta_fname='tdata_R.minor', curiosity_policy=curiosity_policy, max_nodes=max_nodes, regularization=regularization, verbosity=verbosity, map_type=map_type, ablation=ablation, calculate_size=calculate_size, latex_out=latex_out)
    structure(list(rs=rs, rulenames=rulenames, featurenames=featurenames, mat_feature_rule=mat), class="corels")
}

predict.corels <- function(object, tdata, ...) {
    mat_data_feature <- get_data_feature_mat(tdata, object$featurenames)
    mat_data_rules <- mat_data_feature %*% object$mat_feature_rule
    mat_data_rules <- t(t(mat_data_rules)>=c(colSums(object$mat_feature_rule)))+0
    nrules <- ncol(object$mat_feature_rule)
    nsamples <- nrow(tdata)
    mat_idx <- matrix(0, nrow = nsamples, ncol = nrules)

    for (i in 1:length(object$rs$rule)) {
        mat_idx[, object$rs$rule[i]] = i
    }
    mat_satisfy <- mat_data_rules * mat_idx

    # find the earliest rule that captures the data
    mat_caps <- as.matrix(apply(mat_satisfy, 1, function(x) ifelse(!identical(x[x>0], numeric(0)), min(x[x>0]), NaN) ))
    mat_caps[is.na(mat_caps)] = length(object$rs$pred)
    mat_pred <- as.integer(object$rs$pred[mat_caps])
    list(mat_pred)
}

# S3 methods.
# print the model in an interpretable way (if ... then ...)
print.corels <- show.corels <- function(x, useS4 = FALSE, ...) {
    cat(sprintf("\nOPTIMAL RULE LIST\n"))
    for (i in 1:length(x$rs$pred)) {
        if (i==1)
            cat(sprintf("if      %s (rule %d) then predict %d\n", x$rulenames[x$rs$rule[i]], x$rs$rule[i], x$rs$pred[i]))
        else if (i==length(x$rs$pred))
            cat(sprintf("else  (default rule)  then predict %d\n", x$rs$pred[length(x$rs$pred)]))
        else
            cat(sprintf("else if %s (rule %d) then predict %d\n", x$rulenames[x$rs$rule[i]], x$rs$rule[i], x$rs$pred[i]))
    }
}

# This function gets the data-by-feature matrix, given the data and all the feature names
get_data_feature_mat <- function(data, featurenames) {
    mat_data_feature <- matrix(0, nrow=nrow(data), ncol=length(featurenames))
    for (i in 1:length(featurenames)) {
        feature <- featurenames[i]
        conds <- strsplit(feature, '=')
        mat_data_feature[which(data[,conds[[1]][1]]==conds[[1]][2]), i] <- 1
    }
    mat_data_feature
}
