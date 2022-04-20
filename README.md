# flowai operator

##### Description

`flowai` operator performs quality control on flowcytometry data.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the observations (Use 'Time' on top.) 
`y-axis`| measurement value


Output relations|.
---|---
`QCvector`         | numeric, values above 10000 are considered FAIL, under 10000 PASS
`flag`         | character, quality flag `pass` or `fail`

Input parameters|.
---|---
`Detailed` | logical, FALSE will return only the QC_flag. TRUE will also return which event failed what part of the QC, separating the QC flag into flow rate/signal acquisition/dynamic range flags. 
`second_fractionFR`| numeric, fraction of a second that is used to split the time channel in order to recreate the flow rate. Default is 0.1.
`alphaFR` | numeric, level of statistical significance used to accept anomalies detected by the ESD method. Default is 0.01. This setting affects the flow rate QC.
`decompFR` | logical, indicating whether the flow rate should be decomposed in the trend and cyclical components. Default is TRUE and the ESD outlier detection will be executed on the trend component penalized by the magnitude of the cyclical component. If FALSE, the ESD outlier detection will be executed on the original flow rate.
`outlier_binsFS` | logical, indicating whether outlier bins (not events) have to be removed before the changepoint detection of the signal acquisition check. Default is FALSE.
`pen_valueFS` | numeric, value of the penalty for the changepoint detection algorithm. The higher the penalty value the less strict is the detection of the anomalies. Default is 500. This setting affects the signal acquisition QC.
`max_cptFS` | numeric, maximum number of changepoints that can be detected for each channel. Default is 3. This setting affects the signal acquisition QC.
`sideFM` | character string, indicating whether the dynamic range check has to be executed on both limits, the upper limit or the lower limit. Use one of the options: "both", "upper", "lower". Default is "both". 
`neg_valuesFM` | character string, indicating the method to use for the removal of the anomalies from the lower limit of the dynamic range. Use 1 to remove negative outliers or use 2 to truncate the negative values to the cut-off indicated in the FCS file. Default is 1.


##### Details

This operator is able to perform an automatic quality control on FCS data acquired using flow cytometry instruments. By evaluating three different properties: 

* flow rate
* signal acquisition
* dynamic range

The quality control enables the detection and removal of anomalies. The operator returns a QCvector value for each cell, values below 10000 are given to cells who have passed the QC, above 10000 for those who did not pass the QC criteria.

This operator wraps the `flowAI::flow_auto_qc` function from the [flowAI R package](http://bioconductor.org/packages/release/bioc/html/flowAI.html).

##### See Also

[flowsom operator](https://github.com/tercen/flowsom_operator)


