train_from_file.corels <- function(data_fname,label_fname,meta_fname=NULL,curiosity_policy=2,max_nodes=10000, regularization=0.01, verbosity=list("progress"), map_type=2, ablation=0, calculate_size=FALSE, latex_out=FALSE) {
    params_list <- list(curiosity_policy, max_nodes, regularization, verbosity, map_type, ablation, calculate_size, latex_out)
    rs <- .Call('corels_train', PACKAGE='corels', params_list, data_fname, label_fname, meta_fname)
    rs
}

train.corels <- function(tdata,pos_sign="1", neg_sign="0",rule_minlen=1,rule_maxlen=1,minsupport_pos=0.10,minsupport_neg=0.10,curiosity_policy=2,max_nodes=10000, regularization=0.01, verbosity=list("progress"), map_type=2, ablation=0, calculate_size=FALSE, latex_out=FALSE) {
    params_list <- list(curiosity_policy, max_nodes, regularization, verbosity, map_type, ablation, calculate_size, latex_out)
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

    write.table(as.matrix(t(mat_data_rules)), file='tdata_R.out', sep=' ', row.names=rulenames, col.names=FALSE, quote=FALSE)
    label <- t(cbind((tdata$label==neg_sign) +0, (tdata$label==pos_sign) +0))
    write.table(as.matrix(label), file='tdata_R.label', sep=' ', row.names=c("{label=0}", "{label=1}"), col.names=FALSE, quote=FALSE)

    # TODO support for minority points bound

    rs <-.Call('corels_train', PACKAGE='corels', params_list, 'tdata_R.out', 'tdata_R.label', NULL)

    structure(list(rs=rs, rulenames=rulenames, featurenames=featurenames, mat_feature_rule=mat), class="corels")
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
