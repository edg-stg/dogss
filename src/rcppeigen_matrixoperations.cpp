#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd CppEigenInv(Eigen::MatrixXd X) {
  return(X.inverse());
}

// [[Rcpp::export]]
Eigen::MatrixXd CppEigenWoodbury(Eigen::VectorXd V2, Eigen::MatrixXd X, Eigen::VectorXd Sigma2_0) {
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> Lambda=V2.asDiagonal();
  Eigen::MatrixXd XLambda=X*Lambda;
  //Eigen::MatrixXd XLambda=X*V2.asDiagonal();
  //Eigen::MatrixXd Lambda=V2.asDiagonal();
  //Eigen::MatrixXd Inv=;
  Eigen::MatrixXd Vminus=Lambda; //Eigen::MatrixXd(Lambda);
  Vminus.triangularView<Eigen::Upper>() = XLambda.transpose()*((Eigen::MatrixXd(Sigma2_0.asDiagonal())+XLambda*X.transpose()).inverse()).selfadjointView<Eigen::Upper>()*XLambda;
  return(Eigen::MatrixXd(Lambda)-Vminus); //Eigen::MatrixXd(Lambda)
}

// [[Rcpp::export]]
Eigen::VectorXd CppEigenVtimesNu(Eigen::MatrixXd V, Eigen::VectorXd nu) {
  return(V.selfadjointView<Eigen::Upper>()*nu);
}