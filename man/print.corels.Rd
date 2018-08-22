\name{print.corels}
\alias{show.corels}
\alias{print.corels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
INTERPRETABLE VERSION OF A CORELS MODEL
}
\description{
This function prints an corels object. It is a method for the generic function print of class "corels".
}
\usage{
# S3 method for class 'corels'
# This complies with the form of the standard generic method print
\method{print}{corels}(x, useS4=FALSE, \dots)
\method{show}{corels}(x, useS4=FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A CORELS model returned from \code{\link{corels}} function }
  \item{useS4}{ An argument used to match showDefault function. Fixed as FALSE. }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\details{
	 This function is a method for the generic function print for class "corels". It can be invoked by calling print for an object of the appropriate class, or directly by calling print.corels regardless of the class of the object.
}

\examples{
# Let us use the titactoe dataset
data(compas)
for (name in names(compas)) {compas[name] <- as.factor(compas[,name])}

# Train on two-thirds of the data
b = round(2*nrow(compas)/3, digit=0)
data_train <- compas[1:b, ]
# Test on the remaining one third of the data
data_test <- compas[(b+1):nrow(compas), ]
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

