<?php 

/*Evaluation of a linear calibration in analytical chemistry. 
Copyright December 2015, Jose M. Veiga del Baño. 
------------------------------------ 
Recopilation based in the following documents: 
- Handbook of Chemometrics and Qualimetrics. Part A. 
- ISO 8466-2.
- ISO 11843-2.
Distribution probability recopilation of: 
Distribution.php of John Pezullo 
--------------------------------- 
Released under same terms as PHP. 
*/ 
require_once "functions/SchemA.php";
//error_reporting(0);

class secondordercalibration {
  
  //linear calibration var
  var $n; 
  var $X = array(); //array for mass concentration or ratio mass concentration/internal standard (independent variable)
  var $Y = array(); //array for response or ratio response/internal standard (dependent variable)
  var $SumXX;
  var $SumX2Y;
  var $SumXY;
  var $SumX4;
  var $SumX3;
  var $SumX2;
  var $SumX;
  var $SumY;
  var $Den;
  var $SOc;
  var $SOb;
  var $SOa;
  var $Slope;  
  var $PredictedY   = array();
  var $Error        = array();
  var $SquaredError = array();
  var $RSDErrFR     = array(); // Array for response factor
  var $TotalError;  
  var $SumError;
  var $SumSquaredError;  // Residual sum of squares according to ICH and FDA guides
  var $ErrorVariance;
  var $StdErr;          // Standard deviation of the y values in the calibration according to ISO 8466-2
  var $RStdErr;         // Relative standard deviation according to ISO 8466-2
  var $RStdErr100;
  var $Infpoint;        // Test for minima or maxima acording to ISO 8466-2
  var $SumYY;
  var $RSquared;
  var $AlphaTVal;  // T Value for given alpha setting
  var $DeltaTVal;  // T Value for given alpha setting
  var $CCalfa;  // Decision limit according to ISO 11843-2
  var $CCbeta;  //Detection capability according to ISO 11843-2

  function secondordercalibration($X, $Y, $ConfidenceInterval) {

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
    //second-order parameters   
    $this->XMean                = $this->getMean($this->X);
    $this->YMean                = $this->getMean($this->Y);
    $this->SumXX                = $this->getSumXX();
    $this->SumX2Y               = $this->getSumX2Y();
    $this->SumXY                = $this->getSumXY();
    $this->SumX4                = $this->getSumX4();
    $this->SumX3                = $this->getSumX3();
    $this->SumX2                = $this->getSumX2();
    $this->SumX                 = $this->getSumX();
    $this->SumY                 = $this->getSumY();
    $this->Den                  = $this->getDen();
    $this->SOc                  = $this->getSOc(); //Second order c parameter.
    $this->SOb                  = $this->getSOb(); //Second order b parameter.
    $this->SOa                  = $this->getSOa(); //Second order a parameter.
    $this->Slope                = $this->getSlope();
    $this->PredictedY           = $this->getPredictedY();
    $this->Error                = $this->getError();
    $this->SquaredError         = $this->getSquaredError();
    $this->SumError             = $this->getSumError();   
    $this->SumSquaredError      = $this->getSumSquaredError();    
    $this->ErrorVariance        = $this->getErrorVariance();    
    $this->StdErr               = $this->getStdErr();
    $this->RStdErr              = $this->getRStdErr();
    $this->RStdErr100           = $this->getRStdErr100();
    $this->Infpoint             = $this->getInfpoint();
    $this->SumYY                = $this->getSumYY();
    $this->RSquared             = $this->getRSquared();
    $this->CCalfa               = $this->getCCalfa(); 
    $this->CCbeta               = $this->getCCbeta();
    $this->ConfInt              = $ConfidenceInterval;    
    $this->Alpha                = (100 - $this->ConfInt) / 100;
    $this->Delta                = 2*(100 - $this->ConfInt) / 100; //One tail aproximation according to ISO 11843-2
    //Probability parameters 
    $dist = new SchemA;          
    $this->AlphaTVal            = $dist->getInverseStudentT($this->Alpha, $this->X);
    $this->DeltaTVal            = $dist->getInverseStudentT($this->Delta, $this->X);  

    return true;
  }
//***** SECOND ORDER CALIBRATION **** 

//Intermediate calculation for second order calibration
  function getMean($data) {  
    $mean = 0.0;
    $sum  = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $sum += $data[$i];
    }
    $mean  = $sum/$this->n;   
    return $mean;
  }
//Intermediate calculation for second order calibration
  function getSumXX(){  
    $SumXX = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumXX += ($this->X[$i] - $this->XMean) * ($this->X[$i] - $this->XMean);
    }   
    return $SumXX;
  }
  //Intermediate calculation for second order calibration
  function getSumXY(){  
    $SumXY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumXY += ($this->X[$i]) * ($this->Y[$i]);
    }   
    return $SumXY;
  }
  //Intermediate calculation for second order calibration
  function getSumX2Y(){  
    $SumX2Y = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumX2Y += ($this->X[$i]) * ($this->X[$i]) * ($this->Y[$i]);
    }   
    return $SumX2Y;
  }
  //Intermediate calculation for second order calibration
  function getSumX4(){  
    $SumX4 = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumX4 += ($this->X[$i]) * ($this->X[$i])*($this->X[$i])*($this->X[$i]);
    }   
    return $SumX4;
  }
  //Intermediate calculation for second order calibration
  function getSumX3(){  
    $SumX3 = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumX3 += ($this->X[$i]) * ($this->X[$i])*($this->X[$i]);
    }   
    return $SumX3;
  }
  //Intermediate calculation for second order calibration
  function getSumX2(){  
    $SumX2 = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumX2 += ($this->X[$i]) * ($this->X[$i]);
    }   
    return $SumX2;
  }
//Intermediate calculation for second order calibration
  function getSumX(){  
    $SumX = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumX += ($this->X[$i]);
    }   
    return $SumX;
  }
  //Intermediate calculation for second order calibration
  function getSumY(){  
    $SumY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumY += ($this->Y[$i]);
    }   
    return $SumY;
  }
//Intermediate calculation for second order calibration
  function getDen() {
    $Den = 0.0;
    $Den = ($this->n*$this->SumX2*$this->SumX4)+(2*$this->SumX*$this->SumX2*$this->SumX3)-($this->SumX2*$this->SumX2*$this->SumX2)-($this->SumX*$this->SumX*$this->SumX4)-($this->n*$this->SumX3*$this->SumX3);
    return $Den;
  }
//Function to calulate the "c" parameter in a second order calibration acording to ISO 8466-2 (y=a+bx+cx2)
  function getSOc() {
    $SOc = 0.0;
    $SOc = (($this->n*$this->SumX2*$this->SumX2Y)+($this->SumX*$this->SumX3*$this->SumY)+($this->SumX*$this->SumX2*$this->SumXY)-($this->SumX2*$this->SumX2*$this->SumY)-($this->SumX*$this->SumX*$this->SumX2Y)-($this->n*$this->SumX3*$this->SumXY))/$this->Den;
    return $SOc;
  }
//Function to calulate the "b" parameter in a second order calibration acording to ISO 8466-2 (y=a+bx+cx2)
  function getSOb() {
    $SOb = 0.0;
    $SOb = (($this->n*$this->SumX4*$this->SumXY)+($this->SumX*$this->SumX2*$this->SumX2Y)+($this->SumX2*$this->SumX3*$this->SumY)-($this->SumX2*$this->SumX2*$this->SumXY)-($this->SumX*$this->SumX4*$this->SumY)-($this->n*$this->SumX3*$this->SumX2Y))/$this->Den;
    return $SOb;
  }
//Function to calulate parameter the "a" parameter in a second order calibration acording to ISO 8466-2 (y=a+bx+cx2)
  function getSOa() {
    $SOa = 0.0;
    $SOa = (($this->SumX2*$this->SumX4*$this->SumY)+($this->SumX2*$this->SumX3*$this->SumXY)+($this->SumX*$this->SumX3*$this->SumX2Y)-($this->SumX2*$this->SumX2*$this->SumX2Y)-($this->SumX*$this->SumX4*$this->SumXY)-($this->SumX3*$this->SumX3*$this->SumY))/$this->Den;
    return $SOa;
  }
//Function to calculate slope in a second order calibration acording to ISO 5466-2
  function getSlope() {
    $Slope = 0.0;
    $Slope = $this->SOb+($this->XMean*2*$this->SOc);
    return $Slope;
  }
//Intermediate calculation in a second order calibration
  function getPredictedY(){       
    for ($i = 0; $i < $this->n; $i++) {
      $PredictedY[$i] = $this->SOa + ($this->SOb * $this->X[$i])+($this->SOc * ($this->X[$i]* $this->X[$i]));
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
//Intermediate calculation in a second order calibration
  function getSquaredError() {          
    $SquaredError = array();
    for ($i = 0; $i < $this->n; $i++) {
      $SquaredError[$i] = pow(($this->Y[$i] - $this->PredictedY[$i]), 2);
    }   
    return $SquaredError;
  }
//Intermediate calculation in a second order calibration
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
//Intermediate calculation in a second order calibration
  function getErrorVariance() {   
    $ErrorVariance = 0.0;       
    $ErrorVariance = $this->SumSquaredError / ($this->n - 3);   
    return $ErrorVariance;
  }
  //Function to calculate standard deviation of the y values in the calibration according to  ISO 8466-2
  function getStdErr() {   
    $StdErr = 0.0;       
    $StdErr = sqrt($this->ErrorVariance);   
    return $StdErr;
  }
   //Function to calculate relative standard deviation according to  ISO 8466-2
  function getRStdErr() {   
    $RStdErr = 0.0;       
    $RStdErr = (($this->StdErr)/($this->Slope));   
    return $RStdErr;
  }
//Intermediate calculation in a second order calibration
  function getRStdErr100() {   
    $RStdErr100 = 0.0;       
    $RStdErr100 = 100*(($this->RStdErr)/($this->XMean));   
    return $RStdErr100;
  }
//Function to calculate minima or maxima according to  ISO 8466-2
function getInfpoint() {   
    $Infpoint = 0.0;       
    $Infpoint = -(1/2)*($this->SOb)/($this->SOc);   
    return $Infpoint;
  }
//Intermediate calculation in a second order calibration
 function getSumYY(){  
    $SumYY = 0.0;     
    for ($i = 0; $i < $this->n; $i++) {
      $SumYY += ($this->Y[$i] - $this->YMean) * ($this->Y[$i] - $this->YMean);
    }   
    return $SumYY;
  }
//Function to calculate determination coefficient
function getRSquared() {   
    $RSquared = 0.0;       
    $RSquared = 1-(($this->SumError)/($this->SumYY));   
    return $RSquared;
  }
//Function to calculate the decision limit according to ISO 11843-2.
   function getCCalfa() {    
    $CCalfa = 0.0;        
    $CCalfa= ($this->AlphaTVal)* ($this->StdErr/$this->Slope) *sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX));    
    return $CCalfa;
  }
//Function to calculate the detection capability according to ISO 11843-2.
function getCCbeta() {    
    $CCbeta = 0.0;        
    $CCbeta= (2*$this->AlphaTVal)* ($this->StdErr/$this->Slope) *sqrt(1+(1/$this->n)+($this->XMean*$this->XMean)/($this->SumXX));    
    return $CCbeta;
  }
}
