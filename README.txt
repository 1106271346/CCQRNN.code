

Main_function.R is an example to implement the DAU method, IMP method and proposed CCQRNN method.


Before running the R code, please install the package "qrnn", "nleqslv", "quantreg" from https://cran.r-project.org/web/packages/available_packages_by_name.html.


Then run all codes in the Main_function.R, the simulation results for setting (a) are obtained.


Here is the introduction of our R code, to help you know our R code:
  cqr.fit.DA() is used to implemtent the DAU method.
  impute() is the important function to implemtent the IMP method. 
  DArq() is designed to implement the CQRNN method without considering the cure fraction.
  cqr.fit.DAc() is the key function for implementing the proposed CCQRNN method.



Results description:
  hat.gamma:  the estimated cure rate parameters for IMP, DAU, and the CCQRNN method, based on the training set.
  hat.beta:  the estimated parameters associated with survival response for the IMP and DAU methods, based on the training set.
  NN_eta:  the predicted cure rate of the CCQRNN method based on the testing set.
  DA_eta:  the predicted cure rate of the DAU method based on the testing set.
  IMP_eta:  the predicted cure rate of the IMP method based on the testing set.
  MSE_hat:  the predicted MSE of the four methods based on the testing set.
  MAE_hat:  the predicted MAE of the four methods based on the testing set.


Repeat running the main code 100 times to obtain the final simulation results.