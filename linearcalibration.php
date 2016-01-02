<?php 

/*Evaluation of a linear calibration in analytical chemistry. 
Copyright December 2015, Jose M. Veiga del Baño. 
------------------------------------ 
Recopilation based in the following documents: 
- Handbook of Chemometrics and Qualimetrics. Part A. 
- ICH guides. 
- FDA guides.
- CIPAC guides. 
- ISO 8466-1.
- ISO 11095.
- ISO 11843-2.
Distribution probability recopilation of: 
Distribution.php of John Pezullo 
--------------------------------- 
Released under same terms as PHP. 
*/ 
require_once "functions/SchemA.php";
//error_reporting(0);

class linearcalibration {
  
  //linear calibration var
  var $n; 
  var $X = array(); //array for mass concentration or ratio mass concentration/internal standard (independent variable)
  var $Y = array(); //array for response or ratio response/internal standard (dependent variable)
  var $ConfInt;  
  var $Alpha;
  var $XMean;
  var $YMean;
  var $SumXX;
  var $SumXY;
  var $SumYY; 
  var $Slope;
  var $YInt;  
  var $PredictedY   = array();
  var $Error        = array();
  var $SquaredError = array();
  var $RSDErrFR     = array(); // Array for response factor 
  var $SumError;
  var $SumSquaredError;  // Residual sum of squares according to ICH and FDA guides
  var $ErrorVariance;
  var $StdErr;          // Standard deviation of the y values in the calibration according to ISO 8466-1
  var $RStdErr;         // Relative standard deviation according to ISO 8466-1
  var $RStdErr100;
  var $SlopeStdErr;  
  var $YIntStdErr;    
  var $YIntTVal;   // T value for Y Intercept
  var $R;
  var $RSquared;    
  var $DF;         // Degrees of Freedom
  var $SlopeProb;  // Probability of Slope Estimate
  var $YIntProb;   // Probability of Y Intercept Estimate
  var $AlphaTVal;  // T Value for given alpha setting
  var $DeltaTVal;  // T Value for given delta setting
  var $ConfIntOfSlope;
  var $CCalfa;  // Decision limit according to ISO 11843-2
  var $CCbeta;  //Detection capability according to ISO 11843-2
  //weighted linear calibration var
  var $W = array(); //array for Weights in the case of weighted linear calibration.
  var $SumW;
  var $SumWXY;
  var $SumWX;
  var $SumWXX;
  var $SumWYY;
  var $SumWX2;
  var $SumWY2;
  var $SumWY;
  var $SlopeW;
  var $YIntW;
  var $RW;
  var $RSquaredW; 
  var $PredictedYW;
  var $ErrorW;
  var $SquaredErrorW;
  var $SumErrorW;
  var $SumSquaredErrorW;
  var $ErrorVarianceW;
  var $StdErrW;
  var $RStdErrW;
  var $RStdErrW100;
  var $CCalfaW;  // Decision limit according to ISO 11843-2
  var $CCbetaW;  //Detection capability according to ISO 11843-2

  function linearcalibration($X, $Y, $ConfidenceInterval,$W) {

    $numX = count($X);
    $numY = count($Y);

    if ($numX != $numY) {
      die("Error: Size of X and Y vectors must be the same.");

    } 
    if ($numX <= 1) { 
      die("Error: Size of input array must be at least 2.");      
    }
    //General parameters
    $this->n                    = $numX;
    $this->X                    = $X;
    $this->Y                    = $Y;
    $this->W                    = $W;
    //Confidence interval
    $this->ConfInt              = $ConfidenceInterval;    
    $this->Alpha                = (100 - $this->ConfInt) / 100;
    $this->Delta                = (2*(100 - $this->ConfInt)) / 100; //One tail aproximation according to ISO 11843-2
    //Linear parameters   
    $this->XMean                = $this->getMean($this->X);
    $this->YMean                = $this->getMean($this->Y);
    $this->SumXX                = $this->getSumXX();
    $this->SumYY                = $this->getSumYY();
    $this->SumXY                = $this->getSumXY();    
    $this->Slope                = $this->getSlope();
    $this->YInt                 = $this->getYInt();
    $this->PredictedY           = $this->getPredictedY();
    $this->Error                = $this->getError();
    $this->SquaredError         = $this->getSquaredError();
    $this->SumError             = $this->getSumError();   
    $this->SumSquaredError      = $this->getSumSquaredError();    
    $this->ErrorVariance        = $this->getErrorVariance();    
    $this->StdErr               = $this->getStdErr();
    $this->RStdErr              = $this->getRStdErr();
    $this->RStdErr100           = $this->getRStdErr100();
    $this->RSDErrFR             = $this->getRSDErrFR();  
    $this->SlopeStdErr          = $this->getSlopeStdErr();         
    $this->YIntStdErr           = $this->getYIntStdErr();             
    $this->SlopeTVal            = $this->getSlopeTVal();            
    $this->YIntTVal             = $this->getYIntTVal();                
    $this->R                    = $this->getR();                
    $this->RSquared             = $this->getRSquared();  
    $this->DF                   = $this->getDF();
    $this->CCalfa               = $this->getCCalfa(); 
    $this->CCbeta               = $this->getCCbeta(); 
    $this->ConfIntOfSlope       = $this->getConfIntOfSlope();
    //Weight parameters
    $this->SumW                 = $this->getSumW();
    $this->SumWX                = $this->getSumWX();
    $this->SumWXX               = $this->getSumWXX();
    $this->SumWYY               = $this->getSumWYY();
    $this->SumWX2               = $this->getSumWX2();
    $this->SumWY2               = $this->getSumWY2();
    $this->SumWY                = $this->getSumWY();
    $this->SumWXY               = $this->getSumWXY();
    $this->SlopeW               = $this->getSlopeW();
    $this->YIntW                = $this->getYIntW();
    $this->RW                   = $this->getRW();                
    $this->RSquaredW            = $this->getRSquaredW();
    $this->PredictedYW          = $this->getPredictedYW();
    $this->ErrorW               = $this->getErrorW();
    $this->SquaredErrorW        = $this->getSquaredErrorW();
    $this->SumErrorW            = $this->getSumErrorW();   
    $this->SumSquaredErrorW     = $this->getSumSquaredErrorW();    
    $this->ErrorVarianceW       = $this->getErrorVarianceW();    
    $this->StdErrW              = $this->getStdErrW();
    $this->RStdErrW             = $this->getRStdErrW();
    $this->RStdErrW100          = $this->getRStdErrW100();
    $this->CCalfaW              = $this->getCCalfaW(); 
    $this->CCbetaW              = $this->getCCbetaW(); 
    //Probability parameters 
    $dist = new SchemA;  
    $this->SlopeProb            = $dist->getStudentT($this->SlopeTVal, $this->DF);                        
    $this->YIntProb             = $dist->getStudentT($this->YIntTVal, $this->DF);         
    $this->AlphaTVal            = $dist->getInverseStudentT($this->Alpha, $this->X);
    $this->DeltaTVal            = $dist->getInverseStudentT($this->Delta, $this->X);                                  
    
    return true;
  }
//***** LINEAR CALIBRATION **** 

//Intermediate calculation for linear calibration
  function getMean($data) {  
    $mean = 0.0;
    $sum  = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $sum += $data[$i];
    }
    $mean  = $sum/$this->n;   
    return $mean;
  }
//Intermediate calculation for linear calibration
  function getSumXX(){  
    $SumXX = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumXX += ($this->X[$i] - $this->XMean) * ($this->X[$i] - $this->XMean);
    }   
    return $SumXX;
  }
//Intermediate calculation for linear calibration
  function getSumYY(){  
    $SumYY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumYY += ($this->Y[$i] - $this->YMean) * ($this->Y[$i] - $this->YMean);
    }   
    return $SumYY;
  }
//Intermediate calculation for linear calibration
  function getSumXY(){  
    $SumXY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumXY += ($this->X[$i] - $this->XMean) * ($this->Y[$i] - $this->YMean);
    }
    return $SumXY;
  }
//Function to calculate slope in a linear calibration by least squares (ignore origin and no weight)
  function getSlope() {
    $Slope = 0.0;
    $Slope = $this->SumXY / $this->SumXX;
    return $Slope;
  }
//Function to calculate intersection in a linear calibration by least squares (ignore origin and no weight)
  function getYInt() {
    $YInt = 0.0;
    $YInt = $this->YMean - ($this->Slope * $this->XMean);
    return $YInt;
  }
//Intermediate calculation for linear calibration
  function getPredictedY(){       
    for ($i = 0; $i < $this->n; $i++) {
      $PredictedY[$i] = $this->YInt + ($this->Slope * $this->X[$i]);
    }   
    return $PredictedY;
  }
//Function to calculate the array for residual.
  function getError() {          
    $Error = array();
    for ($i = 0; $i < $this->n; $i++) {
      $Error[$i] = $this->Y[$i] - $this->PredictedY[$i];
    }   
    return $Error;
  }
//Intermediate calculation for linear calibration
  function getSquaredError() {          
    $SquaredError = array();
    for ($i = 0; $i < $this->n; $i++) {
      $SquaredError[$i] = pow(($this->Y[$i] - $this->PredictedY[$i]), 2);
    }   
    return $SquaredError;
  }
//Intermediate calculation for linear calibration
  function getSumError() {   //
    $SumError = 0.0;       
    for ($i = 0; $i < $this->n; $i++) {
      $SumError += $this->Error[$i];
    }   
    return $SumError;
  }
  //Function to calculate residual sum of squares according to ICH and FDA guides.
  function getSumSquaredError() {   
    $SumSquaredError = 0.0;       
    for ($i = 0; $i < $this->n; $i++) {
      $SumSquaredError += $this->SquaredError[$i];
    }   
    return $SumSquaredError;
  }
//Intermediate calculation for linear calibration
  function getErrorVariance() {   
    $ErrorVariance = 0.0;       
    $ErrorVariance = $this->SumSquaredError / ($this->n - 2);   
    return $ErrorVariance;
  }
  //Function to calculate standard deviation of the y values in the calibration according to  ISO 8466-1
  function getStdErr() {   
    $StdErr = 0.0;       
    $StdErr = sqrt($this->ErrorVariance);   
    return $StdErr;
  }
  //Function to calculate relative standard deviation according to  ISO 8466-1
  function getRStdErr() {   
    $RStdErr = 0.0;       
    $RStdErr = ($this->StdErr)/($this->Slope);   
    return $RStdErr;
  }
  //Function to calculate relative standard deviation according to  ISO 8466-1
  function getRStdErr100() {   
    $RStdErr100 = 0.0;       
    $RStdErr100 = 100*($this->RStdErr)/($this->XMean);   
    return $RStdErr100;
  }
//Function to calculate array with the response factor (Y versus X). 
  function getRSDErrFR() {         
    for ($i = 0; $i < $this->n; $i++) {
      $RSDErrFR [$i]= ($this->Y[$i])/($this->X[$i]);
    }
    return $RSDErrFR;
  }
//Function to calculate error for the slope according to  Handbook of Chemometric.
  function getSlopeStdErr() {   
    $SlopeStdErr = 0.0;       
    $SlopeStdErr = $this->StdErr / sqrt($this->SumXX);   
    return $SlopeStdErr;
  }
//Function to calculate error for the intercept according to Handbook of Chemometric.
  function getYIntStdErr() {  
    $YIntStdErr = 0.0;       
    $YIntStdErr = $this->StdErr * sqrt(1 / $this->n + pow($this->XMean, 2) / $this->SumXX);   
    return $YIntStdErr;
  }
//Intermediate calculation for linear calibration
  function getSlopeTVal() {   
    $SlopeTVal = 0.0;       
    $SlopeTVal = $this->Slope / $this->SlopeStdErr;   
    return $SlopeTVal;
  }
//Intermediate calculation
  function getYIntTVal() {   
    $YIntTVal = 0.0;       
    $YIntTVal = $this->YInt / $this->YIntStdErr;  
    return $YIntTVal;
  }
//Function to calculate correlation coefficient
  function getR() {   
    $R = 0.0;       
    $R = $this->SumXY / sqrt($this->SumXX * $this->SumYY);
    return $R;
  }
//Function to calculate determination coefficient
  function getRSquared() {  
    $RSquared = 0.0;       
    $RSquared = $this->R * $this->R;   
    return $RSquared;
  }
//Function to calculate freedom degrees
  function getDF() {    
    $DF = 0.0;       
    $DF = $this->n - 2;
    return $DF;
  }
//Function to calculate determination coefficient
  function getConfIntOfSlope() {    
    $ConfIntOfSlope = 0.0;        
    $ConfIntOfSlope = $this->AlphaTVal * $this->SlopeStdErr ;        
    return $ConfIntOfSlope;
  }
//Function to calculate the critical value according to ISO 11843-2 in a linear calibration
   function getCCalfa() {    
    $CCalfa = 0.0;        
    $CCalfa= ($this->YInt)+ (($this->AlphaTVal)* ($this->StdErr) *sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX)));    
    return $CCalfa;
  }
//Function to calculate the minimun detectable value according to ISO 11843-2 in a linear calibration
function getCCbeta() {    
    $CCbeta = 0.0;        
    $CCbeta= 2*($this->DeltaTVal)* ($this->StdErr/$this->Slope)*sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX));       
    return $CCbeta;
  }
  
//***** WEIGHTED LINEAR CALIBRATION  ****

//Intermediate calculation for weighted linear calibration
  function getSumW(){  
    $SumW = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumW += $this->W[$i];
    }   
    return $SumW;
  }
  //Intermediate calculation for weighted linear calibration
  function getSumWXY(){  
    $SumWXY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWXY += ($this->W[$i]*$this->X[$i]*$this->Y[$i]);
    }
    return $SumWXY;
  }
   //Intermediate calculation for weighted linear calibration
  function getSumWX(){  
    $SumWX = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWX += ($this->W[$i]*$this->X[$i]);
    }
    return $SumWX;
  }
//Intermediate calculation for weighted linear calibration
function getSumWY(){  
    $SumWY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWY += ($this->W[$i]*$this->Y[$i]);
    }
    return $SumWY;
  }
//Intermediate calculation for weighted linear calibration
function getSumWXX(){  
    $SumWXX = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWXX += ($this->W[$i]*$this->X[$i]*$this->X[$i]);
    }
    return $SumWXX;
  }
  function getSumWYY(){  
    $SumWYY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWYY += ($this->W[$i]*$this->Y[$i]*$this->Y[$i]);
    }
    return $SumWYY;
  }
//Intermediate calculation for weighted linear calibration
  function getSumWX2(){  
    $SumWX2 = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWX2 += ($this->W[$i]*$this->X[$i]);
    }
    return $SumWX2*$SumWX2;
  }
//Intermediate calculation for weighted linear calibration
    function getSumWY2(){  
    $SumWY2 = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumWY2 += ($this->W[$i]*$this->Y[$i]);
    }
    return $SumWY2*$SumWY2;
  }
//Function to calculate slope in a weighted linear calibration by least squares (ignore origin)
  function getSlopeW() {
    $SlopeW = 0.0;
    $SlopeW = (($this->SumW*$this->SumWXY)-($this->SumWX*$this->SumWY))/(($this->SumW*$this->SumWXX)-($this->SumWX2));
    return $SlopeW;
  }
  //Function to calculate intersection in a weighted linear calibration by least squares (ignore origin)
  function getYIntW() {
    $YIntW = 0.0;
    $YIntW = (($this->SumWXX*$this->SumWY)-($this->SumWX*$this->SumWXY))/(($this->SumW*$this->SumWXX)-($this->SumWX2));;
    return $YIntW;
  }
  //Function to calculate correlation coefficient in a weighted linear calibration
  function getRW() {   
    $RW = 0.0;       
    $RW = (($this->SumW*$this->SumWXY)-($this->SumWX*$this->SumWY))/(sqrt(($this->SumW*$this->SumWXX)-($this->SumWX2))*sqrt(($this->SumW*$this->SumWYY)-($this->SumWY2)));
    return $RW;
  }
//Function to calculate determination coefficient in a weighted linear calibration
  function getRSquaredW() {  
    $RSquaredW = 0.0;       
    $RSquaredW = $this->RW * $this->RW;   
    return $RSquaredW;
  }
 //Intermediate calculation for linear calibration in a weighted linear calibration
  function getPredictedYW(){       
    for ($i = 0; $i < $this->n; $i++) {
      $PredictedYW[$i] = $this->YIntW + ($this->SlopeW * $this->X[$i]);
    }   
    return $PredictedYW;
  }
//Function to calculate the array for residual in a weighted linear calibration
  function getErrorW() {          
    $ErrorW = array();
    for ($i = 0; $i < $this->n; $i++) {
      $ErrorW[$i] = $this->Y[$i] - $this->PredictedYW[$i];
    }   
    return $ErrorW;
  }
//Intermediate calculation for linear calibration in a weighted linear calibration
  function getSquaredErrorW() {          
    $SquaredErrorW = array();
    for ($i = 0; $i < $this->n; $i++) {
      $SquaredErrorW[$i] = pow(($this->Y[$i] - $this->PredictedYW[$i]), 2);
    }   
    return $SquaredErrorW;
  }
//Intermediate calculation for linear calibration in a weighted linear calibration
  function getSumErrorW() {   //
    $SumErrorW = 0.0;       
    for ($i = 0; $i < $this->n; $i++) {
      $SumErrorW += $this->ErrorW[$i];
    }   
    return $SumErrorW;
  }
  //Function to calculate residual sum of squares according to ICH and FDA guides in a weighted linear calibration
  function getSumSquaredErrorW() {   
    $SumSquaredErrorW = 0.0;       
    for ($i = 0; $i < $this->n; $i++) {
      $SumSquaredErrorW += $this->SquaredErrorW[$i];
    }   
    return $SumSquaredErrorW;
  }
//Intermediate calculation for linear calibration in a weighted linear calibration
  function getErrorVarianceW() {   
    $ErrorVarianceW = 0.0;       
    $ErrorVarianceW = $this->SumSquaredErrorW / ($this->n - 2);   
    return $ErrorVarianceW;
  }
  //Function to calculate standard deviation of the y values in the calibration according to  ISO 8466-1 in a weighted linear calibration
  function getStdErrW() {   
    $StdErrW = 0.0;       
    $StdErrW = sqrt($this->ErrorVarianceW);   
    return $StdErrW;
  }
   //Function to calculate relative standard deviation according to  ISO 8466-1 in a weighted linear calibration
  function getRStdErrW() {   
    $RStdErrW = 0.0;       
    $RStdErrW = ($this->StdErrW)/($this->SlopeW);   
    return $RStdErrW;
  }
  //Function to calculate relative standard deviation according to  ISO 8466-1 in a weighted linear calibration
  function getRStdErrW100() {   
    $RStdErrW100 = 0.0;       
    $RStdErrW100 = 100*($this->RStdErrW)/($this->XMean);   
    return $RStdErrW100;
  }
  //Function to calculate the critical value according to ISO 11843-2 in a weighted linear calibration
   function getCCalfaW() {    
    $CCalfaW = 0.0;        
    $CCalfaW= ($this->YIntW)+ (($this->AlphaTVal)* ($this->StdErrW) *sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX)));      
    return $CCalfaW;
  }
//Function to calculate the minimun detectable value according to ISO 11843-2 in a weighted linear calibration
function getCCbetaW() {    
    $CCbetaW = 0.0;        
    $CCbetaW= 2*($this->DeltaTVal)* ($this->StdErrW/$this->SlopeW)*sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX));   
    return $CCbetaW;
  }
}
