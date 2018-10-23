\name{corels-package}
\alias{corels-package}
\docType{package}
\title{
  CORELS: LEARNING CERTIFIABLY OPTIMAL RULE LISTS
}
\description{
Fit a \code{corels} model.
Learn from the data and create a decision rule list of binary predictions in the format of:

  if (condition1) then predict ...

  else if (condition2) then predict ...

  else if (condition3) then predict ...

  ...

  else (default rule) then predict ...

( See the examples below )
}
\details{
This package contains three functions: \code{\link{train.corels}}, \code{\link{print.corels}}, and \code{\link{predict.corels}}
}
\author{
Aditya Saligrama, Elaine Angelino, Nicholas Larus-Stone, Daniel Alabi, Vassilios Kaxiras, Cynthia Rudin, Margo Seltzer
}
\references{
  Elaine Angelino, Nicholas Larus-Stone, Daniel Alabi, Margo Seltzer and Cynthia Rudin (2017)
  \emph{Learning Certifiably Optimal Rule Lists for Categorical Data.}
  Journal of Machine Learning Research, 2018.

  Nicholas Larus-Stone, Elaine Angelino, Daniel Alabi, Margo Seltzer, Vassilios Kaxiras, Aditya Saligrama and Cynthia Rudin (2018)
  \emph{Systems Optimizations for Learning Certifiably Optimal Rule Lists.}
  SysML Conference, 2018.
}
\keyword{ package }
\seealso{
\code{\link{train.corels}}, \code{\link{print.corels}}, and \code{\link{predict.corels}}
}
\examples{
# Let us use the tictactoe dataset
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
# Restrict max nodes to 10000 so it doesn't take too much time
# For most cases, increase it to 10^8 or greater
corels_model <- train.corels(data_train,rule_maxlen=2,
  minsupport_pos=0.005,minsupport_neg=0.005,max_nodes=10000)

print(corels_model)

# Make predictions on the test set
yhat <- predict(corels_model, data_test)
# yhat will be a list of binary predictions for the test data.
}
