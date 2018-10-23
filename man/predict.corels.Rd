\name{predict.corels}
\alias{predict}
\alias{predict.corels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
PREDICT OUTCOMES FOR EACH OBSERVATION IN TEST DATA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a list of binary predictions.
}
\usage{
\method{predict}{corels}(object, tdata, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{corels model returned from the \code{\link{train.corels}} function.
%%     ~~Describe \code{model} here~~
}
  \item{tdata}{test data
%%     ~~Describe \code{tdata} here~~
}
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\value{
Returns a list containing binary predictions for each observation in the test data according to the rule list model.
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
# Restrict max nodes to 10000 so it doesn't take too much time
# For most cases, increase it to 10^8 or greater
corels_model <- train.corels(data_train,rule_maxlen=2,
  minsupport_pos=0.005,minsupport_neg=0.005,max_nodes=10000)

print(corels_model)

# Make predictions on the test set
  yhat <- predict(corels_model, data_test)
# yhat will be a list of binary predictions for the test data.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
