<?php 

/*Common Statistical Analysis for Chemical laboratories. 
Copyright December 2015, Jose M. Veiga del Baño. 
------------------------------------ 
Web recopilation and creation of the following documents: 
- Handbook of Chemometrics and Qualimetrics. Part A. 
- NORDTEST 569. 
- NORDTEST 537. 
- ISO 13528. 
Distribution probability recopilation of: 
Distribution.php of John Pezullo 
--------------------------------- 
Released under same terms as PHP. 
*/ 
class SchemA { 
  
  var $n; 
  var $mean; 
  var $percentile; 
  var $probability; 
  var $df; 
  var $reference; 
  var $X = array(); 
  var $Y = array(); 

//Function to calculate arithmetic average 
function average($X){ 
    if (!count($X)) return 0; 
    $sum = 0; 
    for ($i = 0; $i < count($X); $i++) 
    { 
        $sum += $X[$i]; 
    } 

    return $sum / count($X); 
} 
// Function to calculate standard deviation (n-1) 
function sd($X) { 
    if (count($X)<2) return 0; 
    $aver=$this->average($X); 
 $sum = 0; 
    for ($i = 0; $i < count($X); $i++) 
    { 
        $sum += ($X[$i]-$aver)*($X[$i]-$aver); 
    } 
    return sqrt($sum / (count($X)-1)); 
} 
//Function to calculate percentiles 
function mypercentile($X,$percentile){ 
    if( 0 < $percentile && $percentile < 1 ) { 
        $p = $percentile; 
    }else if( 1 < $percentile && $percentile <= 100 ) { 
        $p = $percentile * .01; 
    }else { 
        return ""; 
    } 
    $count = count($X); 
    $allindex = ($count-1)*$p; 
    $intvalindex = intval($allindex); 
    $floatval = $allindex - $intvalindex; 
    sort($X); 
    if(!is_float($floatval)){ 
        $result = $X[$intvalindex]; 
    }else { 
        if($count > $intvalindex+1) 
            $result = $floatval*($X[$intvalindex+1] - $X[$intvalindex]) + $X[$intvalindex]; 
        else 
            $result = $X[$intvalindex]; 
    } 
    return $result; 
} 
//Function to calculate Scaled median absolute deviation (MAD). 
function mad($X){ 
    if (!count($X)) return 0; 
    $median = $this->mypercentile($X,50); 
        if (count($X)<10){ 
       $sum = 0; 
       for ($i = 0; $i < count($X); $i++) 
            { 
            $sum += abs($X[$i]-$median); 
            } 

        $result =(1/(0.798*count($X)))*$sum; 
        return $result; 
        }else{ 
        for ($i = 0; $i < count($X); $i++) 
            { 
            $medi[]= abs($X[$i]-$median); 
            } 
        $med=$this->mypercentile($medi,50); 
        $result =1.483*$med; 
        return $result; 
    } 
} 
//Function to calculate relative standard deviation (RSD). Coefficient of variation is 100*RSD. 
function rsd($X){ 
    $result =($this->sd($X))/($this->average($X)); 
    return $result; 
} 
//Function to calculate relative standard deviation (RSD) with robust estimators. Coefficient of variation is 100*RSD. 
function rsdrob($X){ 
    $result =($this->mad($X))/($this->mypercentile($X,50)); 
    return $result; 
} 
//Function to calculate percentage of root mean square (RMS) for a reference value (CRM, recovery,...). 
function rms($X,$reference){ 
    $sesgo=0; 
    for($i=0;$i<count($X);$i++){ 
        $sesgo+=($X[$i]-$reference)*($X[$i]-$reference); 
    } 
    $raiz=sqrt($sesgo / count($X)); 
return $raiz; 
} 
 //Function for F and t distributions 
  function doCommonMath($q, $i, $j, $b) { 
    
    $zz = 1; 
    $z = $zz; 
    $k = $i; 
    
    
    while($k <= $j) { 
      $zz = $zz * $q * $k / ($k - $b); 
      $z = $z + $zz; 
      $k = $k + 2; 
    } 
    return $z; 
  } 
 //Function to calculate the value of the distribution F 
 function getFisherF($f, $n1, $n2) { 
    
    $x = $n2 / ($n1 * $f + $n2); 
        
    if(($n1%2)==0) { 
      return $this->doCommonMath(1-$x, $n2, $n1+$n2-4, $n2-2) * pow($x, $n2/2); 
    } 
    if(($n2%2)==0){ 
      return 1 - $this->doCommonMath($x, $n1, $n1+$n2-4, $n1-2) * pow(1-$x, $n1/2); 
    } 
    $th = atan(sqrt($n1 * $f / $n2)); 
    $a = $th / (pi() / 2); 
    $sth = sin($th); 
    $cth = cos($th); 
    if($n2 > 1) { 
      $a = $a + $sth * $cth * $this->doCommonMath($cth*$cth, 2, $n2-3, -1) / (pi()/2); 
    } 
    if($n1==1) { 
      return 1 - $a; 
    } 
    $c = 4 * $this->doCommonMath($sth*$sth, $n2+1, $n1+$n2-4, $n2-2)* $sth * pow($cth,$n2) / pi(); 
    if($n2==1) { 
      return 1 - $a + $c/2; 
    } 
    $k=2; 
    while($k<=($n2-1)/2) { 
      $c = $c * $k/($k-.5); 
      $k=$k+1; 
    } 
    return 1-$a+$c; 
  } 

  //Function to calculate the inverse value of the distribution F 
  function getInverseFisherF($probability, $X, $Y) { 
     if( 0 < $probability && $probability < 1 ) { 
        $p = $probability; 
    }else { 
        $p = (100-$probability) * 0.01; 
    } 
    $n1=count($X)-1; 
    $n2=count($Y)-1; 
    $v = 0.5; 
    $dv = 0.5; 
    $f = 0.0; 
   
    while($dv > 1e-10) { 
      
      $f = (1 / $v) - 1; 
      $dv = $dv / 2; 

      if($this->getFisherF($f, $n1, $n2) > $p) { 
        $v = $v - $dv; 
      } else { 
        $v = $v + $dv; 
      } 
    } 
    return $f; 
  } 
  //Function to calculate the value of the distribution t-student 
  function getStudentT($t, $df) { 

    $t = abs($t); 
    $w = $t / sqrt($df); 
    $th = atan($w); 
    
    if ($df == 1) { 
      return 1 - $th / (pi() / 2); 
    } 
  
    $sth = sin($th); 
    $cth = cos($th); 
  
    if( ($df % 2) ==1 ) { 
      return 1 - ($th + $sth * $cth * $this->doCommonMath($cth * $cth, 2, $df - 3, -1)) / (pi()/2); 
    } else { 
      return 1 - $sth * $this->doCommonMath($cth * $cth, 1, $df - 3, -1); 
    } 
  
  } 
  //Function to calculate the inverse value of the distribution t-student 
  function getInverseStudentT($probability, $X) { 
    if( 0 < $probability && $probability < 1 ) { 
        $p = $probability; 
    }else { 
        $p = (100-$probability) * 0.01; 
    } 
    $df=count($X)-1; 
    $v = 0.5; 
    $dv = 0.5; 
    $t = 0; 
    
    while($dv > 1e-6) { 
      $t = (1 / $v) - 1; 
      $dv = $dv / 2; 
      if ( $this->getStudentT($t, $df) > $p) { 
        $v = $v - $dv; 
      } else { 
        $v = $v + $dv; 
      } 
    } 
    return $t; 
  } 
//Function to calculate the value of F for two populations 
 function getFcalc($X, $Y) { 
   $variancex=$this->sd($X)*$this->sd($X); 
   $variancey=$this->sd($Y)*$this->sd($Y); 
   $maxim=max($variancex,$variancey); 
   $minim=min($variancex,$variancey); 
   $fcalc=$maxim/$minim; 
   return $fcalc; 
  } 

//Function to calculate the value of t for a reference value (CRM, recovery..) 
 function getTcalc($X, $reference) { 
   $avervalues=$this->average($X); 
   $sdvalues=$this->sd($X); 
   $n=count($X); 
   $tcalc=(abs($avervalues-$reference)*sqrt($n))/$sdvalues; 
   return $tcalc; 
  } 
  /*Estimation of measurement relative uncertainty based in the following assumptions: 
  - The data $X are in reproducibility conditions (RSD). The laboratory not analyzed in repeatability conditions. 
  - The outlier datas are identified and removed. 
  - The reference value is the theorical value or 100 for recovery data. 
  - RMS is the only one estimation for the bias of the measurement. 
  - Coberture factor is the t-student value for a 95% of probability. 
  */ 
  function uncertainty($X,$reference) { 
   $median=$this->mypercentile($X,50); 
   $mad=$this->mad($X); 
   $outliermax=$median+3*$mad; 
   $outliermin=$median-3*$mad; 
   for ($i = 0; $i < count($X); $i++) 
    { 
        if( $X[$i]> $outliermin & $X[$i]< $outliermax ){ 
          $Z[] = $X[$i]; 
        } 
    } 
   $reprod=$this->rsd($Z); 
   $bias=$this->rms($Z,$reference); 
   $df=count($Z)-1; 
   $coberture=$this->getInverseStudentT(95,$Z); 
   $uncert=100*($coberture*sqrt($reprod*$reprod+$bias*$bias)); 
   return $uncert; 
  } 
//Function to formats a number ($value) to a specified number of significant figures. 
  function sigFig($value, $sigFigs) { 
  $exponent = floor(log10(abs($value))+1); 
  $significand = round(($value 
    / pow(10, $exponent)) 
    * pow(10, $sigFigs)) 
    / pow(10, $sigFigs); 
  return $significand * pow(10, $exponent); 
} 

} 
