\name{train.corels}
\alias{train.corels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
TRAIN THE CORELS MODEL WITH THE GIVEN DATA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
TRAIN THE CORELS MODEL WITH THE GIVEN DATA
}
\usage{
train.corels(tdata, pos_sign="1", neg_sign="0", rule_minlen=1, rule_maxlen=1, minsupport_pos=0.10, minsupport_neg=0.10, curiosity_policy=2, max_nodes=10000, regularization = 0.01, verbosity="progress", map_type=1, ablation=0, calculate_size=0, latex_out=1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tdata}{a dataframe, with a "label" column specifying the correct labels for each observation.}

  \item{pos_sign}{the sign for the positive labels in the "label" column.}

  \item{neg_sign}{the sign for the negative labels in the "label" column.}

  \item{rule_minlen}{the minimum number of cardinality for rules to be mined from the dataframe.}

  \item{rule_maxlen}{the maximum number of cardinality for rules to be mined from the dataframe.}

  \item{minsupport_pos}{a number between 0 and 1, for the minimum percentage support for the positive observations.}

  \item{minsupport_neg}{a number between 0 and 1, for the minimum percentage support for the negative observations.}

  \item{curiosity_policy}{1 to use the custom curiosity ordering, 2 to use lower bound ordering, 3 to use objective ordering, 4 to use pure DFS, or 5 to use pure BFS.}

  \item{max_nodes}{maximum number of nodes in the CORELS queue. Use to limit runtime at the cost of a possibly non-optimal rule list.}

  \item{verbosity}{a comma-separated string comprised of "progress", "log", "rule", "label", "samples", "silent". If the silent option is used, it must be the only option provided.}

  \item{map_type}{0 to use no map, 1 to use the prefix permutation map, or 2 to use the captured symmetry map.}

  \item{ablation}{0 to use all optimizations, 1 to exclude the minimum support bounds, 2 to exclude the lookahead bound,}

  \item{calculate_size}{calculate upper bound on remaining search space.}
}
\value{
Return a list of :
	\item{rs}{a ruleset which contains the rule indices and their binary predictions for the best rule list by training corels with the given data and parameters.}

	\item{rulenames}{a list of all the rule names mined with \code{arules}.}

	\item{featurenames}{a list of all the feature names.}

	\item{mat_feature_rule}{a binary matrix representing which features are included in which rules.}
}
\references{
  Elaine Angelino, Nicholas Larus-Stone, Daniel Alabi, Margo Seltzer and Cynthia Rudin (2017)
  \emph{Learning Certifiably Optimal Rule Lists for Categorical Data.}
  Journal of Machine Learning Research, 2018.

  Nicholas Larus-Stone, Elaine Angelino, Daniel Alabi, Margo Seltzer, Vassilios Kaxiras, Aditya Saligrama and Cynthia Rudin (2018)
  \emph{Systems Optimizations for Learning Certifiably Optimal Rule Lists.}
  SysML Conference, 2018.
}
\author{
%%  ~~who you are~~
Aditya Saligrama, Elaine Angelino, Nicholas Larus-Stone, Daniel Alabi, Vassilios Kaxiras, Cynthia Rudin, Margo Seltzer
}
\examples{
# Let us use the titactoe dataset
data(tictactoe)
for (name in names(tictactoe)) {tictactoe[name] <- as.factor(tictactoe[,name])}

# Train on two-thirds of the data
b = round(2*nrow(tictactoe)/3, digit=0)
data_train <- tictactoe[1:b, ]
# Test on the remaining one third of the data
data_test <- tictactoe[(b+1):nrow(tictactoe), ]
# data_train, data_test are dataframes with factor columns
# The class column is "label"

# Run the CORELS algorithm on the training set
corels_model <- train.corels(data_train,rule_maxlen=2,minsupport_pos=0.005,minsupport_neg=0.005,max_nodes=10000000000)
print(corels_model)

# Make predictions on the test set
yhat <- predict(corels_model, data_test)
# yhat will be a list of binary predictions for the test data.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
